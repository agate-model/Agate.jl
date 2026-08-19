using ForwardDiff

using Agate.Compilation:
    GrowthParameterBinding,
    GrowthBiomassContribution,
    GrowthResourceLossContribution,
    realize_process_topology,
    process_contributions,
    group_contributions,
    compile_tendencies
using Agate.Configuration: realize_components
using Agate.Construction: define_tracer_functions
using Agate.Equations: CompiledEquation
using Agate.Factories: default_components, default_community
using Agate.Processes: ModelDefinition, normalize_model

const NIPIZD_GROWTH_TARGET_ORDER = (:N, :P_1, :P_2)
const NIPIZD_GROWTH_ARGS = (
    0.0, 0.0, 0.0, 0.0, 7.0, 1.0, 0.05, 0.08, 0.01, 0.02, 100.0
)

function nipizd_growth_compilation(::Type{T}=Float64) where {T<:Real}
    factory = Agate.Models.NiPiZD.NiPiZDFactory()
    normalized = normalize_model(ModelDefinition(factory))
    layout = realize_components(default_components(factory); scalar_type=T)
    context = Agate.Configuration.parse_community(
        T, default_community(factory); biogeochem_tracers=(:N, :D)
    )
    process = normalized.processes.growth_P
    topology = realize_process_topology(process, layout, context)
    binding = GrowthParameterBinding(normalized, :growth_P)
    contributions = process_contributions(process, topology, binding)
    grouped = group_contributions(contributions; target_order=NIPIZD_GROWTH_TARGET_ORDER)
    return (; binding, topology, contributions, grouped, compiled=compile_tendencies(grouped), context)
end

@testset "Smith/Monod growth contribution compiler" begin
    compilation = nipizd_growth_compilation()
    binding = compilation.binding
    topology = compilation.topology
    contributions = compilation.contributions
    grouped = compilation.grouped
    compiled = compilation.compiled

    @test binding.maximum_rate === :maximum_growth_rate
    @test binding.half_saturation === :nutrient_half_saturation
    @test binding.alpha === :alpha
    @test topology.population_tracers == (:P_1, :P_2)
    @test topology.population_indices == (3, 4)
    @test topology.resource_target === :N
    @test topology.light_driver === :PAR

    @test length(contributions) == 4
    @test count(c -> c isa GrowthBiomassContribution, contributions) == 2
    @test count(c -> c isa GrowthResourceLossContribution, contributions) == 2
    @test map(length, grouped) == (N=2, P_1=1, P_2=1)
    @test keys(compiled) == NIPIZD_GROWTH_TARGET_ORDER
    @test all(equation -> equation isa CompiledEquation, values(compiled))
    @test all(equation -> isbitstype(typeof(equation.f)), values(compiled))
    @test all(
        equation -> all(term -> isbitstype(typeof(term)), equation.f.terms), values(compiled)
    )

    bgc = Agate.Models.NiPiZD.construct(;
        parameters=(
            linear_mortality=(Z_1=0.0, Z_2=0.0, P_1=0.0, P_2=0.0),
            quadratic_mortality=(Z_1=0.0, Z_2=0.0),
            maximum_growth_rate=(P_1=2.2e-5, P_2=1.7e-5),
            nutrient_half_saturation=(P_1=0.21, P_2=0.34),
            alpha=(P_1=2.4e-6, P_2=1.8e-6),
            maximum_predation_rate=(Z_1=0.0, Z_2=0.0),
            detritus_remineralization=0.0,
        ),
    )
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
        target -> generated_bgc(Val(target), NIPIZD_GROWTH_ARGS...),
        NIPIZD_GROWTH_TARGET_ORDER,
    )
    legacy = map(
        target -> bgc(Val(target), NIPIZD_GROWTH_ARGS...),
        NIPIZD_GROWTH_TARGET_ORDER,
    )
    @test all(process_compiler_isapprox.(generated, legacy))
    @test isapprox(sum(generated), 0; atol=10 * eps(sum(abs, generated)))
    @test process_compiler_isapprox(
        @inferred(compiled.P_1(bgc, NIPIZD_GROWTH_ARGS...)), legacy[2]
    )
end

@testset "Smith/Monod growth compiler ForwardDiff" begin
    compiled = nipizd_growth_compilation().compiled
    bgc = Agate.Models.NiPiZD.construct()
    biomass = 0.013

    function p_1_growth(P)
        args = (
            0.0, 0.0, 0.0, 0.0, 7.0, 1.0, 0.05, 0.08, P, 0.02, 100.0
        )
        return compiled.P_1(bgc, args...)
    end

    rate = p_1_growth(biomass)
    derivative = ForwardDiff.derivative(p_1_growth, biomass)
    @test process_compiler_isapprox(derivative, rate / biomass)
end
