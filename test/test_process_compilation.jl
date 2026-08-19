using ForwardDiff
using Enzyme

using Agate.Compilation:
    MortalityParameterBinding,
    MortalityLossContribution,
    MortalityRoutingContribution,
    realize_process_topology,
    process_contributions,
    group_contributions,
    compile_tendencies
using Agate.Configuration: realize_components
using Agate.Construction: define_tracer_functions
using Agate.Equations: CompiledEquation
using Agate.Factories: default_components, default_community
using Agate.Processes: ModelDefinition, normalize_model

const NIPIZD_MORTALITY_IDS = (
    :linear_mortality_P, :linear_mortality_Z, :quadratic_mortality_Z
)
const NIPIZD_TRACER_ORDER = (:N, :D, :Z_1, :Z_2, :P_1, :P_2)
const NIPIZD_MORTALITY_ARGS = (
    0.0, 0.0, 0.0, 0.0, 7.0, 1.0, 0.05, 0.08, 0.01, 0.02, 100.0
)

function nipizd_mortality_compilation(::Type{T}=Float64) where {T<:Real}
    factory = Agate.Models.NiPiZD.NiPiZDFactory()
    normalized = normalize_model(ModelDefinition(factory))
    layout = realize_components(default_components(factory); scalar_type=T)
    context = Agate.Configuration.parse_community(
        T, default_community(factory); biogeochem_tracers=(:N, :D)
    )
    contributions = ()
    for id in NIPIZD_MORTALITY_IDS
        process = getproperty(normalized.processes, id)
        topology = realize_process_topology(process, layout, context)
        binding = MortalityParameterBinding(normalized, id)
        contributions = (
            contributions..., process_contributions(process, topology, binding)...
        )
    end
    grouped = group_contributions(contributions; target_order=NIPIZD_TRACER_ORDER)
    return (; contributions, grouped, compiled=compile_tendencies(grouped), context)
end

@testset "Mortality process contribution compiler" begin
    compilation = nipizd_mortality_compilation()
    contributions = compilation.contributions
    grouped = compilation.grouped
    compiled = compilation.compiled

    normalized = normalize_model(
        ModelDefinition(Agate.Models.NiPiZD.NiPiZDFactory())
    )
    @test MortalityParameterBinding(normalized, :linear_mortality_P).rate ===
        :linear_mortality
    @test MortalityParameterBinding(normalized, :quadratic_mortality_Z).rate ===
        :quadratic_mortality
    @test MortalityParameterBinding(
        normalized, :linear_mortality_Z
    ).routing_fraction === :mortality_export_fraction

    @test length(contributions) == 18
    @test count(contribution -> contribution isa MortalityLossContribution, contributions) == 6
    @test count(contribution -> contribution isa MortalityRoutingContribution, contributions) == 12
    loss_indices(process) = Tuple(
        contribution.population_index for contribution in contributions if
        contribution isa MortalityLossContribution && contribution.process === process
    )
    @test loss_indices(:linear_mortality_P) == (3, 4)
    @test loss_indices(:linear_mortality_Z) == (1, 2)
    @test map(length, grouped) == (N=6, D=6, Z_1=2, Z_2=2, P_1=1, P_2=1)
    @test keys(compiled) == NIPIZD_TRACER_ORDER
    @test all(equation -> equation isa CompiledEquation, values(compiled))
    @test all(equation -> isbitstype(typeof(equation.f)), values(compiled))
    @test all(
        equation -> all(term -> isbitstype(typeof(term)), equation.f.terms), values(compiled)
    )

    bgc = Agate.Models.NiPiZD.construct(;
        parameters=(
            linear_mortality=(Z_1=1.1e-6, Z_2=1.3e-6, P_1=0.7e-6, P_2=0.9e-6),
            quadratic_mortality=(Z_1=1.7e-6, Z_2=2.1e-6),
            maximum_growth_rate=(P_1=0.0, P_2=0.0),
            maximum_predation_rate=(Z_1=0.0, Z_2=0.0),
            detritus_remineralization=0.0,
        ),
    )
    args = NIPIZD_MORTALITY_ARGS
    tracer_index = Agate.Runtime.build_tracer_index(
        compilation.context,
        NIPIZD_TRACER_ORDER,
        (:PAR,);
        n_biogeochem_tracers=2,
    )
    factory = define_tracer_functions(
        bgc.parameters, compiled; auxiliary_fields=(:PAR,), tracer_index
    )
    generated_bgc = factory(bgc.parameters)
    @test !any(type -> type === Any, fieldtypes(typeof(generated_bgc.tracer_functions)))

    for target in NIPIZD_TRACER_ORDER
        generated = generated_bgc(Val(target), args...)
        legacy = bgc(Val(target), args...)
        @test process_compiler_isapprox(generated, legacy)
    end
    @test process_compiler_isapprox(
        @inferred(compiled.P_1(bgc, args...)), bgc(Val(:P_1), args...)
    )
end

@testset "Mortality compiler ForwardDiff" begin
    compiled = nipizd_mortality_compilation().compiled
    biomass = 0.013

    function p_1_mortality(rate)
        T = typeof(rate)
        bgc = Agate.Models.NiPiZD.construct(;
            scalar_type=T, parameters=(linear_mortality=(P_1=rate,),)
        )
        args = (
            zero(T), zero(T), zero(T), zero(T), T(7), T(1), T(0.05), T(0.05),
            T(biomass), T(0.02), T(100),
        )
        return compiled.P_1(bgc, args...)
    end

    rate = 0.9e-6
    derivative = ForwardDiff.derivative(p_1_mortality, rate)
    @test process_compiler_isapprox(derivative, -biomass)
end

@testset "Mortality compiler Enzyme" begin
    compiled = nipizd_mortality_compilation().compiled
    base_bgc = Agate.Models.NiPiZD.construct(;
        parameters=(linear_mortality=(P_1=0.9e-6,), mortality_export_fraction=0.2)
    )
    active = Agate.Runtime.active_parameters(
        base_bgc; linear_mortality=(:P_1,), mortality_export_fraction=true
    )
    args = NIPIZD_MORTALITY_ARGS

    function diagnostic(parameters)
        bgc = Agate.Runtime.parameterized(base_bgc, parameters; active_parameters=active)
        return compiled.P_1(bgc, args...) + 0.25 * compiled.D(bgc, args...) +
               0.5 * compiled.N(bgc, args...)
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
    @testset "Mortality compiler GPU execution" begin
        @eval using CUDA
        @eval using OceanBioME: Biogeochemistry
        @eval using Oceananigans: RectilinearGrid, NonhydrostaticModel, set!, time_step!
        @eval using Oceananigans.Architectures: GPU
        @eval using Oceananigans.Grids: Periodic, Flat, Bounded

        @test CUDA.functional()
        if CUDA.functional()
            compilation = nipizd_mortality_compilation(Float32)
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
                NIPIZD_TRACER_ORDER,
                ();
                n_biogeochem_tracers=2,
            )
            factory = define_tracer_functions(
                base_bgc.parameters,
                compilation.compiled;
                auxiliary_fields=(),
                tracer_index,
            )
            mortality_bgc = factory(base_bgc.parameters)
            bgc_model = Biogeochemistry(mortality_bgc)
            model = NonhydrostaticModel(; grid, biogeochemistry=bgc_model)
            set!(
                model;
                N=7f0,
                D=1f0,
                Z_1=0.05f0,
                Z_2=0.08f0,
                P_1=0.01f0,
                P_2=0.02f0,
            )
            time_step!(model, 60f0)
            @test model.clock.iteration == 1
        end
    end
end
