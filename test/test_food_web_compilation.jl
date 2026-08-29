using ForwardDiff
using Oceananigans.Biogeochemistry:
    required_biogeochemical_auxiliary_fields, required_biogeochemical_tracers

using Agate.Configuration: Plankton, Pool
using Agate.Construction: construct
using Agate.Introspection: interaction_matrix
using Agate.Parameters: Parameter, NoDefault
using Agate.Processes:
    ModelDefinition, Growth, Light, NutrientResponse, Temperature, Consumption, Smith, Monod,
    Q10, HeterotrophicConsumption, PreferentialGrazing, participants

function food_web_definition()
    components = (
        N=Pool(:nitrogen),
        D=Pool(:nitrogen),
        POM=Pool(:nitrogen),
        P=Plankton(; states=:nitrogen, reference_state=:nitrogen, size_structure=[1.0]),
        B=Plankton(; states=:nitrogen, reference_state=:nitrogen, size_structure=[0.8]),
        M=Plankton(; states=:nitrogen, reference_state=:nitrogen, size_structure=[2.0]),
        Z=Plankton(; states=:nitrogen, reference_state=:nitrogen, size_structure=[10.0]),
    )
    temperature = Temperature(
        Q10(); bindings=(q10=:temperature_q10, reference_temperature=:reference_temperature)
    )
    processes = (
        growth_autotrophs=Growth(;
            plankton=(:P, :M),
            bindings=(maximum_rate=:maximum_growth_rate,),
            factors=(
                temperature=temperature,
                nutrients=NutrientResponse(
                    Monod(); resource=:N,
                    bindings=(half_saturation=:nutrient_half_saturation,)
                ),
                light=Light(Smith(); driver=:PAR),
            ),
        ),
        consume_POM=Consumption(
            HeterotrophicConsumption();
            consumers=:B,
            resources=:POM,
            bindings=(
                maximum_rate=:maximum_consumption_rate,
                half_saturation=:pom_half_saturation,
                assimilation=:bacterial_assimilation,
            ),
            factors=(temperature=temperature,),
            unassimilated_products=:D,
        ),
        grazing_living=Consumption(
            PreferentialGrazing();
            consumers=(:M, :Z),
            resources=(:P, :B),
            bindings=(
                maximum_rate=:maximum_predation_rate,
                half_saturation=:holling_half_saturation,
                palatability=:living_palatability_matrix,
                assimilation=:living_assimilation_matrix,
            ),
            unassimilated_products=:D,
        ),
    )
    no_default() = Parameter(NoDefault())
    parameters = (
        maximum_growth_rate=no_default(),
        alpha=no_default(),
        nutrient_half_saturation=no_default(),
        temperature_q10=no_default(),
        reference_temperature=no_default(),
        maximum_consumption_rate=no_default(),
        pom_half_saturation=no_default(),
        bacterial_assimilation=no_default(),
        maximum_predation_rate=no_default(),
        holling_half_saturation=no_default(),
        living_palatability_matrix=no_default(),
        living_assimilation_matrix=no_default(),
    )
    return ModelDefinition(; components, processes, parameters)
end

function food_web_parameter_overrides(::Type{T}=Float64) where {T<:Real}
    return (
        maximum_growth_rate=T[2e-5, 1.4e-5],
        alpha=T[2e-6, 1.6e-6],
        nutrient_half_saturation=T[0.2, 0.3],
        temperature_q10=T(2),
        reference_temperature=T(20),
        maximum_consumption_rate=T[1.5e-5],
        pom_half_saturation=T[0.15],
        bacterial_assimilation=reshape(T[0.65], 1, 1),
        maximum_predation_rate=T[6e-5, 9e-5],
        holling_half_saturation=T[0.12, 0.18],
        living_palatability_matrix=T[0.6 0.8; 0.7 0.9],
        living_assimilation_matrix=T[0.4 0.5; 0.35 0.45],
    )
end

function food_web_args(bgc, state::NamedTuple; PAR=0.0, temperature=20.0)
    tracers = required_biogeochemical_tracers(bgc)
    tracer_values = Tuple(
        hasproperty(state, tracer) ? getproperty(state, tracer) : 0.0 for tracer in tracers
    )
    auxiliary_values = Tuple(
        auxiliary === :PAR ? PAR :
        auxiliary === :temperature ? temperature :
        error("unknown test auxiliary field :$auxiliary")
        for auxiliary in required_biogeochemical_auxiliary_fields(bgc)
    )
    return (0.0, 0.0, 0.0, 0.0, tracer_values..., auxiliary_values...)
end

@testset "POM, bacteria, mixotrophy, and reusable factors" begin
    definition = food_web_definition()
    bgc = construct(definition; parameter_overrides=food_web_parameter_overrides())

    @test participants(definition.processes.consume_POM) == (
        consumer=(:B,), resource=(:POM,)
    )
    @test participants(definition.processes.grazing_living).resource == (:P, :B)
    @test :POM ∉ participants(definition.processes.grazing_living).resource
    @test :M ∈ participants(definition.processes.growth_autotrophs).plankton
    @test :M ∈ participants(definition.processes.grazing_living).consumer
    @test required_biogeochemical_auxiliary_fields(bgc) == (:PAR, :temperature)

    state = (
        N=5.0, D=0.1, POM=0.5,
        P_1=0.05, B_1=0.03, M_1=0.02, Z_1=0.04,
    )
    args = food_web_args(bgc, state; PAR=100.0, temperature=25.0)
    tendencies = values(model_tendencies(bgc, args))
    @test isapprox(sum(tendencies), 0; atol=10 * eps(sum(abs, tendencies)))

    consumption_state = (POM=0.5, B_1=0.03)
    consumption20 = bgc(
        Val(:POM), food_web_args(bgc, consumption_state; temperature=20.0)...
    )
    consumption30 = bgc(
        Val(:POM), food_web_args(bgc, consumption_state; temperature=30.0)...
    )
    @test process_compiler_isapprox(consumption30, 2 * consumption20)

    growth_state = (N=5.0, P_1=0.05)
    growth20 = bgc(
        Val(:P_1), food_web_args(bgc, growth_state; PAR=100.0, temperature=20.0)...
    )
    growth30 = bgc(
        Val(:P_1), food_web_args(bgc, growth_state; PAR=100.0, temperature=30.0)...
    )
    @test process_compiler_isapprox(growth30, 2 * growth20)

    derivative = ForwardDiff.derivative(0.5) do pom
        dynamic_state = (POM=pom, B_1=0.03)
        bgc(Val(:POM), food_web_args(bgc, dynamic_state; temperature=25.0)...)
    end
    @test isfinite(derivative)
    @test derivative < 0
end

@testset "Multi-resource consumer storage axes" begin
    components = (
        N=Pool(:nitrogen),
        POM_1=Pool(:nitrogen),
        POM_2=Pool(:nitrogen),
        POM_3=Pool(:nitrogen),
        X=Plankton(; states=:nitrogen, reference_state=:nitrogen, size_structure=[0.4]),
        B=Plankton(; states=:nitrogen, reference_state=:nitrogen, size_structure=[0.8]),
    )
    processes = (
        consume_POM=Consumption(
            HeterotrophicConsumption();
            consumers=:B,
            resources=(:POM_1, :POM_2, :POM_3),
            bindings=(
                maximum_rate=:maximum_consumption_rate,
                half_saturation=:pom_half_saturation,
                assimilation=:bacterial_assimilation,
            ),
            unassimilated_products=:N,
        ),
    )
    parameters = (
        maximum_consumption_rate=Parameter(NoDefault()),
        pom_half_saturation=Parameter(NoDefault()),
        bacterial_assimilation=Parameter(NoDefault()),
    )
    definition = ModelDefinition(; components, processes, parameters)
    bgc = construct(
        definition;
        parameter_overrides=(
            maximum_consumption_rate=[2.0],
            pom_half_saturation=[1.0, 3.0, 7.0],
            bacterial_assimilation=reshape([0.2, 0.4, 0.8], 1, 3),
        ),
    )

    @test bgc.parameters.maximum_consumption_rate == [2.0]
    @test bgc.parameters.pom_half_saturation == [1.0, 3.0, 7.0]
    @test bgc.parameters.bacterial_assimilation == reshape([0.2, 0.4, 0.8], 1, 3)
    @test_throws ArgumentError interaction_matrix(bgc, :bacterial_assimilation)

    names = Agate.Introspection.tracer_names(bgc)
    state = (N=0.0, POM_1=1.0, POM_2=1.0, POM_3=1.0, X_1=5.0, B_1=1.0)
    args = (0.0, 0.0, 0.0, 0.0, Tuple(getproperty(state, name) for name in names)...)
    expected = (POM_1=-1.0, POM_2=-0.5, POM_3=-0.25, B_1=0.6, N=1.15, X_1=0.0)

    for (name, value) in pairs(expected)
        @test bgc(Val(name), args...) ≈ value
    end
    @test isapprox(sum(bgc(Val(name), args...) for name in names), 0.0; atol=1e-14)
end
