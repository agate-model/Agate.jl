using ForwardDiff
using Enzyme

using Agate.Compilation:
    RemineralizationParameterBinding,
    RemineralizationSourceLossContribution,
    RemineralizationDestinationGainContribution,
    realize_process_topology,
    process_contributions,
    model_contributions,
    group_contributions,
    compile_model_tendencies
using Agate.Configuration: realize_components
using Agate.Construction: define_tracer_functions
using Agate.Equations: CompiledEquation
using Agate.Factories: default_components, default_community
using Agate.Processes: ModelDefinition, driver_identities, normalize_model

const NIPIZD_PROCESS_TRACER_ORDER = (:N, :D, :Z_1, :Z_2, :P_1, :P_2)
const NIPIZD_PROCESS_ARGS = (
    0.0, 0.0, 0.0, 0.0, 7.0, 1.0, 0.05, 0.08, 0.01, 0.02, 100.0
)

function nipizd_process_compilation(::Type{T}=Float64) where {T<:Real}
    factory = Agate.Models.NiPiZD.NiPiZDFactory()
    normalized = normalize_model(ModelDefinition(factory))
    layout = realize_components(default_components(factory); scalar_type=T)
    context = Agate.Configuration.parse_community(
        T, default_community(factory); biogeochem_tracers=(:N, :D)
    )
    contributions = model_contributions(normalized, layout, context)
    grouped = group_contributions(
        contributions; target_order=NIPIZD_PROCESS_TRACER_ORDER
    )
    compiled = compile_model_tendencies(
        normalized, layout, context; target_order=NIPIZD_PROCESS_TRACER_ORDER
    )
    drivers = driver_identities(normalized)
    return (; normalized, layout, context, contributions, grouped, compiled, drivers)
end

function full_process_bgc()
    return Agate.Models.NiPiZD.construct(;
        parameters=(
            detritus_remineralization=1.4e-6,
            mortality_export_fraction=0.23,
            linear_mortality=(Z_1=1.1e-6, Z_2=1.3e-6, P_1=0.7e-6, P_2=0.9e-6),
            quadratic_mortality=(Z_1=1.7e-6, Z_2=2.1e-6),
            maximum_growth_rate=(P_1=2.2e-5, P_2=1.7e-5),
            nutrient_half_saturation=(P_1=0.21, P_2=0.34),
            alpha=(P_1=2.4e-6, P_2=1.8e-6),
            maximum_predation_rate=(Z_1=1.1e-4, Z_2=0.8e-4),
            holling_half_saturation=(Z_1=0.12, Z_2=0.18),
        ),
        palatability_matrix=[0.7 0.2; 0.4 0.9],
        assimilation_matrix=[0.3 0.35; 0.4 0.45],
    )
end

@testset "Linear remineralization contribution compiler" begin
    compilation = nipizd_process_compilation()
    normalized = compilation.normalized
    named = normalized.processes.remineralization_D
    binding = RemineralizationParameterBinding(normalized, :remineralization_D)
    topology = realize_process_topology(named, compilation.layout, compilation.context)
    contributions = process_contributions(named, topology, binding)

    @test binding.rate === :detritus_remineralization
    @test topology.source_tracers == (:D,)
    @test topology.destination_target === :N
    @test length(contributions) == 2
    @test contributions[1] isa RemineralizationSourceLossContribution
    @test contributions[2] isa RemineralizationDestinationGainContribution
    @test contributions[1].target === :D
    @test contributions[2].target === :N
end

@testset "Complete NiPiZD process contribution compiler" begin
    compilation = nipizd_process_compilation()
    contributions = compilation.contributions
    grouped = compilation.grouped
    compiled = compilation.compiled

    @test compilation.drivers == (:PAR,)
    @test length(contributions) == 36
    @test map(length, grouped) == (N=9, D=11, Z_1=4, Z_2=4, P_1=4, P_2=4)
    @test keys(compiled) == NIPIZD_PROCESS_TRACER_ORDER
    @test all(equation -> equation isa CompiledEquation, values(compiled))
    @test all(equation -> isbitstype(typeof(equation.f)), values(compiled))
    @test all(
        equation -> all(term -> isbitstype(typeof(term)), equation.f.terms), values(compiled)
    )

    bgc = full_process_bgc()
    tracer_index = Agate.Runtime.build_tracer_index(
        compilation.context,
        NIPIZD_PROCESS_TRACER_ORDER,
        compilation.drivers;
        n_biogeochem_tracers=2,
    )
    factory = define_tracer_functions(
        bgc.parameters,
        compiled;
        auxiliary_fields=compilation.drivers,
        tracer_index,
    )
    generated_bgc = factory(bgc.parameters)
    @test !any(type -> type === Any, fieldtypes(typeof(generated_bgc.tracer_functions)))

    generated = map(
        target -> generated_bgc(Val(target), NIPIZD_PROCESS_ARGS...),
        NIPIZD_PROCESS_TRACER_ORDER,
    )
    constructed = map(
        target -> bgc(Val(target), NIPIZD_PROCESS_ARGS...),
        NIPIZD_PROCESS_TRACER_ORDER,
    )
    @test all(process_compiler_isapprox.(generated, constructed))
    @test isapprox(sum(generated), 0; atol=10 * eps(sum(abs, generated)))
    @test process_compiler_isapprox(
        @inferred(compiled.P_1(bgc, NIPIZD_PROCESS_ARGS...)), constructed[5]
    )
end

@testset "Complete NiPiZD compiler ForwardDiff" begin
    compiled = nipizd_process_compilation().compiled
    bgc = full_process_bgc()
    biomass = 0.013

    args(P) = (0.0, 0.0, 0.0, 0.0, 7.0, 1.0, 0.05, 0.08, P, 0.02, 100.0)
    generated_tendency(P) = compiled.P_1(bgc, args(P)...)
    constructed_tendency(P) = bgc(Val(:P_1), args(P)...)

    generated_derivative = ForwardDiff.derivative(generated_tendency, biomass)
    constructed_derivative = ForwardDiff.derivative(constructed_tendency, biomass)
    @test process_compiler_isapprox(generated_derivative, constructed_derivative)
end

@testset "Complete NiPiZD compiler Enzyme" begin
    compiled = nipizd_process_compilation().compiled
    base_bgc = full_process_bgc()
    active = Agate.Runtime.active_parameters(
        base_bgc;
        maximum_growth_rate=(:P_1,),
        linear_mortality=(:P_1,),
        detritus_remineralization=true,
    )
    args = NIPIZD_PROCESS_ARGS

    function diagnostic(parameters)
        bgc = Agate.Runtime.parameterized(base_bgc, parameters; active_parameters=active)
        return compiled.P_1(bgc, args...) + 0.5 * compiled.Z_1(bgc, args...) +
               0.25 * compiled.D(bgc, args...) + 0.125 * compiled.N(bgc, args...)
    end

    p0 = copy(active.values)
    gradient = zeros(length(p0))
    Enzyme.autodiff(
        Enzyme.set_runtime_activity(Enzyme.Reverse),
        Enzyme.Const(diagnostic),
        Enzyme.Active,
        Enzyme.Duplicated(p0, gradient),
    )
    @test all(isfinite, gradient)

    for i in eachindex(p0)
        delta = max(abs(p0[i]), 1.0) * 1e-6
        plus = copy(p0)
        minus = copy(p0)
        plus[i] += delta
        minus[i] -= delta
        finite_difference = (diagnostic(plus) - diagnostic(minus)) / (2delta)
        @test isapprox(gradient[i], finite_difference; rtol=1e-4, atol=1e-10)
    end
end

if lowercase(get(ENV, "AGATE_TEST_CUDA", "0")) in ("1", "true", "yes")
    @testset "Complete NiPiZD compiler GPU execution" begin
        @eval using CUDA
        @eval using Agate.Library.Light: CyclicalPAR, FunctionFieldPAR
        @eval using OceanBioME: Biogeochemistry
        @eval using Oceananigans: RectilinearGrid, NonhydrostaticModel, set!, time_step!
        @eval using Oceananigans.Architectures: GPU
        @eval using Oceananigans.Grids: Periodic, Flat, Bounded

        @test CUDA.functional()
        if CUDA.functional()
            compilation = nipizd_process_compilation(Float32)
            grid = RectilinearGrid(
                GPU(), Float32;
                topology=(Periodic, Flat, Bounded),
                size=(4, 4),
                x=(0f0, 4f0),
                z=(-4f0, 0f0),
            )
            base_bgc = Agate.Models.NiPiZD.construct(; grid)
            tracer_index = Agate.Runtime.build_tracer_index(
                compilation.context,
                NIPIZD_PROCESS_TRACER_ORDER,
                compilation.drivers;
                n_biogeochem_tracers=2,
            )
            factory = define_tracer_functions(
                base_bgc.parameters,
                compilation.compiled;
                auxiliary_fields=compilation.drivers,
                tracer_index,
            )
            compiled_bgc = factory(base_bgc.parameters)
            light_attenuation = FunctionFieldPAR(; grid, PAR_f=CyclicalPAR())
            model = NonhydrostaticModel(;
                grid,
                biogeochemistry=Biogeochemistry(compiled_bgc; light_attenuation),
            )
            set!(
                model;
                N=7f0,
                D=0.01f0,
                Z_1=0.05f0,
                Z_2=0.05f0,
                P_1=0.01f0,
                P_2=0.01f0,
            )
            time_step!(model, 60f0)
            @test model.clock.iteration == 1
        end
    end
end
