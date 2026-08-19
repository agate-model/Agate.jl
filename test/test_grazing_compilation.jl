using ForwardDiff

using Agate.Compilation:
    GrazingParameterBinding,
    GrazingResourceLossContribution,
    GrazingConsumerGainContribution,
    GrazingUnassimilatedContribution,
    realize_process_topology,
    process_contributions,
    group_contributions,
    compile_tendencies
using Agate.Configuration: realize_components
using Agate.Construction: define_tracer_functions
using Agate.Equations: CompiledEquation
using Agate.Factories: default_components, default_community
using Agate.Processes: ModelDefinition, normalize_model

const NIPIZD_GRAZING_TARGET_ORDER = (:D, :Z_1, :Z_2, :P_1, :P_2)
const NIPIZD_GRAZING_ARGS = (
    0.0, 0.0, 0.0, 0.0, 7.0, 1.0, 0.05, 0.08, 0.01, 0.02, 100.0
)
const NIPIZD_TEST_PALATABILITY = [0.7 0.2; 0.4 0.9]
const NIPIZD_TEST_ASSIMILATION = [0.3 0.35; 0.4 0.45]

function nipizd_grazing_compilation(::Type{T}=Float64) where {T<:Real}
    factory = Agate.Models.NiPiZD.NiPiZDFactory()
    normalized = normalize_model(ModelDefinition(factory))
    layout = realize_components(default_components(factory); scalar_type=T)
    context = Agate.Configuration.parse_community(
        T, default_community(factory); biogeochem_tracers=(:N, :D)
    )
    process = normalized.processes.grazing_Z_on_P
    topology = realize_process_topology(process, layout, context)
    binding = GrazingParameterBinding(normalized, :grazing_Z_on_P)
    contributions = process_contributions(process, topology, binding)
    grouped = group_contributions(contributions; target_order=NIPIZD_GRAZING_TARGET_ORDER)
    return (; binding, topology, contributions, grouped, compiled=compile_tendencies(grouped), context)
end

function grazing_only_bgc()
    return Agate.Models.NiPiZD.construct(;
        parameters=(
            linear_mortality=(Z_1=0.0, Z_2=0.0, P_1=0.0, P_2=0.0),
            quadratic_mortality=(Z_1=0.0, Z_2=0.0),
            maximum_growth_rate=(P_1=0.0, P_2=0.0),
            maximum_predation_rate=(Z_1=1.1e-4, Z_2=0.8e-4),
            holling_half_saturation=(Z_1=0.12, Z_2=0.18),
            detritus_remineralization=0.0,
        ),
        palatability_matrix=NIPIZD_TEST_PALATABILITY,
        assimilation_matrix=NIPIZD_TEST_ASSIMILATION,
    )
end

@testset "Preferential grazing contribution compiler" begin
    compilation = nipizd_grazing_compilation()
    binding = compilation.binding
    topology = compilation.topology
    contributions = compilation.contributions
    grouped = compilation.grouped
    compiled = compilation.compiled

    @test binding.maximum_rate === :maximum_predation_rate
    @test binding.half_saturation === :holling_half_saturation
    @test binding.palatability === :palatability_matrix
    @test binding.assimilation === :assimilation_matrix

    @test topology.consumer_tracers == (:Z_1, :Z_2)
    @test topology.consumer_indices == (1, 2)
    @test topology.resource_tracers == (:P_1, :P_2)
    @test topology.resource_indices == (3, 4)
    @test topology.unassimilated_target === :D
    @test compilation.context.consumer_indices == collect(1:4)
    @test compilation.context.prey_indices == collect(1:4)

    @test length(contributions) == 12
    @test count(c -> c isa GrazingResourceLossContribution, contributions) == 4
    @test count(c -> c isa GrazingConsumerGainContribution, contributions) == 4
    @test count(c -> c isa GrazingUnassimilatedContribution, contributions) == 4
    resource_pairs = Tuple(
        (c.consumer_axis, c.resource_axis) for c in contributions if
        c isa GrazingResourceLossContribution
    )
    @test resource_pairs == ((1, 1), (1, 2), (2, 1), (2, 2))
    @test map(length, grouped) == (D=4, Z_1=2, Z_2=2, P_1=2, P_2=2)
    @test keys(compiled) == NIPIZD_GRAZING_TARGET_ORDER
    @test all(equation -> equation isa CompiledEquation, values(compiled))
    @test all(equation -> isbitstype(typeof(equation.f)), values(compiled))
    @test all(
        equation -> all(term -> isbitstype(typeof(term)), equation.f.terms), values(compiled)
    )

    bgc = grazing_only_bgc()
    tracer_index = Agate.Runtime.build_tracer_index(
        compilation.context,
        (:N, :D, :Z_1, :Z_2, :P_1, :P_2),
        (:PAR,);
        n_biogeochem_tracers=2,
    )
    factory = define_tracer_functions(
        bgc.parameters, compiled; auxiliary_fields=(:PAR,), tracer_index
    )
    generated_bgc = factory(bgc.parameters)
    @test !any(type -> type === Any, fieldtypes(typeof(generated_bgc.tracer_functions)))

    generated = map(
        target -> generated_bgc(Val(target), NIPIZD_GRAZING_ARGS...),
        NIPIZD_GRAZING_TARGET_ORDER,
    )
    legacy = map(
        target -> bgc(Val(target), NIPIZD_GRAZING_ARGS...),
        NIPIZD_GRAZING_TARGET_ORDER,
    )
    @test all(process_compiler_isapprox.(generated, legacy))
    @test isapprox(sum(generated), 0; atol=10 * eps(sum(abs, generated)))
    @test process_compiler_isapprox(
        @inferred(compiled.P_1(bgc, NIPIZD_GRAZING_ARGS...)), legacy[4]
    )
end

@testset "Preferential grazing compiler ForwardDiff" begin
    compiled = nipizd_grazing_compilation().compiled
    bgc = grazing_only_bgc()
    biomass = 0.013

    args(P) = (0.0, 0.0, 0.0, 0.0, 7.0, 1.0, 0.05, 0.08, P, 0.02, 100.0)
    generated_loss(P) = compiled.P_1(bgc, args(P)...)
    legacy_loss(P) = bgc(Val(:P_1), args(P)...)

    generated_derivative = ForwardDiff.derivative(generated_loss, biomass)
    legacy_derivative = ForwardDiff.derivative(legacy_loss, biomass)
    @test process_compiler_isapprox(generated_derivative, legacy_derivative)
    @test generated_derivative < 0
end

if lowercase(get(ENV, "AGATE_TEST_CUDA", "0")) in ("1", "true", "yes")
    @testset "Preferential grazing compiler GPU execution" begin
        @eval using CUDA
        @eval using OceanBioME: Biogeochemistry
        @eval using Oceananigans: RectilinearGrid, NonhydrostaticModel, set!, time_step!
        @eval using Oceananigans.Architectures: GPU
        @eval using Oceananigans.Grids: Periodic, Flat, Bounded

        @test CUDA.functional()
        if CUDA.functional()
            compilation = nipizd_grazing_compilation(Float32)
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
                NIPIZD_GRAZING_TARGET_ORDER,
                ();
                n_biogeochem_tracers=1,
            )
            factory = define_tracer_functions(
                base_bgc.parameters,
                compilation.compiled;
                auxiliary_fields=(),
                tracer_index,
            )
            grazing_bgc = factory(base_bgc.parameters)
            model = NonhydrostaticModel(;
                grid, biogeochemistry=Biogeochemistry(grazing_bgc)
            )
            set!(model; D=1f0, Z_1=0.05f0, Z_2=0.08f0, P_1=0.01f0, P_2=0.02f0)
            time_step!(model, 60f0)
            @test model.clock.iteration == 1
        end
    end
end
