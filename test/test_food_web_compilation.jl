using ForwardDiff
using Oceananigans.Biogeochemistry:
    required_biogeochemical_auxiliary_fields, required_biogeochemical_tracers

using Agate.Components: Plankton, Pool
using Agate.Construction: construct
using Agate.Parameters: Parameter, NoDefault
using Agate.Processes:
    ModelDefinition, Growth, Light, NutrientResponse, Temperature, Consumption, Smith, Monod,
    Q10, HeterotrophicConsumption, PreferentialGrazing, participants

function food_web_definition(; grazing=PreferentialGrazing())
    components = (
        N=Pool(:nitrogen),
        D=Pool(:nitrogen),
        POM=Pool(:nitrogen),
        P=Plankton(; states=(nitrogen=:nitrogen,), reference_state=:nitrogen, size_structure=[1.0]),
        B=Plankton(; states=(nitrogen=:nitrogen,), reference_state=:nitrogen, size_structure=[0.8]),
        M=Plankton(; states=(nitrogen=:nitrogen,), reference_state=:nitrogen, size_structure=[2.0]),
        Z=Plankton(; states=(nitrogen=:nitrogen,), reference_state=:nitrogen, size_structure=[10.0]),
    )
    temperature = Temperature(
        Q10(); bindings=(q10=:temperature_q10, reference_temperature=:reference_temperature)
    )
    processes = (
        growth_autotrophs=Growth(;
            plankton=(:P, :M),
            reference_resource=:N,
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
                substrate_preference=:substrate_preference_matrix,
                assimilation=:bacterial_assimilation,
            ),
            factors=(temperature=temperature,),
            unassimilated_products=:D,
        ),
        grazing_living=Consumption(
            grazing;
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
        substrate_preference_matrix=no_default(),
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
        substrate_preference_matrix=reshape(T[1.0], 1, 1),
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
    valid = food_web_parameter_overrides()
    for (name, value, domain, shown) in (
        (:maximum_growth_rate, [NaN, 1.4e-5], :nonnegative, "NaN"),
        (:maximum_growth_rate, [-1.0, 1.4e-5], :nonnegative, "-1.0"),
        (:reference_temperature, Inf, :finite, "Inf"),
        (:temperature_q10, 0.0, :positive, "0.0"),
        (:living_palatability_matrix, [NaN 0.8; 0.7 0.9], :nonnegative, "NaN"),
        (:living_assimilation_matrix, [-0.1 0.5; 0.35 0.45], :unit_interval, "-0.1"),
        (:living_assimilation_matrix, [1.1 0.5; 0.35 0.45], :unit_interval, "1.1"),
    )
        overrides = merge(valid, NamedTuple{(name,)}((value,)))
        message = argument_error_message(() -> construct(definition; parameter_overrides=overrides))
        @test all(occursin(fragment, message) for fragment in
            ("process :", "parameter :$name", "domain :$domain", shown))
    end

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
    direct_consumption20 = -1.5e-5 * 0.03 *
        Agate.Processes.factor_value(Monod(), 0.5, 0.15)
    @test process_compiler_isapprox(consumption20, direct_consumption20)

    growth_state = (N=5.0, P_1=0.05)
    growth20 = bgc(
        Val(:P_1), food_web_args(bgc, growth_state; PAR=100.0, temperature=20.0)...
    )
    growth30 = bgc(
        Val(:P_1), food_web_args(bgc, growth_state; PAR=100.0, temperature=30.0)...
    )
    @test process_compiler_isapprox(growth30, 2 * growth20)
    direct_growth20 = 0.05 * 2e-5 *
        Agate.Processes.factor_value(Q10(), 20.0, 2.0, 20.0) *
        Agate.Processes.factor_value(Monod(), 5.0, 0.2) *
        Agate.Processes.factor_value(Smith(), 100.0, 2e-5, 2e-6)
    @test process_compiler_isapprox(growth20, direct_growth20)

    derivative = ForwardDiff.derivative(0.5) do pom
        dynamic_state = (POM=pom, B_1=0.03)
        bgc(Val(:POM), food_web_args(bgc, dynamic_state; temperature=25.0)...)
    end
    @test isfinite(derivative)
    @test derivative < 0
end

@testset "Multi-resource heterotrophs share capacity across substrates" begin
    components = (
        N=Pool(:nitrogen),
        POM_1=Pool(:nitrogen),
        POM_2=Pool(:nitrogen),
        POM_3=Pool(:nitrogen),
        B=Plankton(; states=(nitrogen=:nitrogen,), reference_state=:nitrogen, size_structure=[0.8]),
    )
    processes = (
        consume_POM=Consumption(
            HeterotrophicConsumption();
            consumers=:B,
            resources=(:POM_1, :POM_2, :POM_3),
            bindings=(
                maximum_rate=:maximum_consumption_rate,
                half_saturation=:pom_half_saturation,
                substrate_preference=:substrate_preference_matrix,
                assimilation=:bacterial_assimilation,
            ),
            unassimilated_products=:N,
        ),
    )
    parameters = (
        maximum_consumption_rate=Parameter(NoDefault()),
        pom_half_saturation=Parameter(NoDefault()),
        substrate_preference_matrix=Parameter(NoDefault()),
        bacterial_assimilation=Parameter(NoDefault()),
    )
    definition = ModelDefinition(; components, processes, parameters)
    base_overrides = (
        maximum_consumption_rate=[2.0],
        pom_half_saturation=[1.0, 3.0, 7.0],
        substrate_preference_matrix=ones(1, 3),
        bacterial_assimilation=reshape([0.2, 0.4, 0.8], 1, 3),
    )
    bgc = construct(definition; parameter_overrides=base_overrides)

    names = Agate.Introspection.tracer_names(bgc)
    state = (N=0.0, POM_1=1.0, POM_2=1.0, POM_3=1.0, B_1=1.0)
    args = (0.0, 0.0, 0.0, 0.0, Tuple(getproperty(state, name) for name in names)...)
    expected = (POM_1=-21 / 26, POM_2=-7 / 26, POM_3=-3 / 26, B_1=47 / 130, N=54 / 65)

    for (name, value) in pairs(expected)
        @test bgc(Val(name), args...) ≈ value
    end
    total_uptake = -sum(bgc(Val(name), args...) for name in (:POM_1, :POM_2, :POM_3))
    @test total_uptake < 2.0
    @test isapprox(sum(bgc(Val(name), args...) for name in names), 0.0; atol=1e-14)

    preferred = construct(
        definition;
        parameter_overrides=merge(
            base_overrides,
            (substrate_preference_matrix=reshape([1.0, 3.0, 7.0], 1, 3),),
        ),
    )
    preferred_losses = Tuple(-preferred(Val(name), args...) for name in (:POM_1, :POM_2, :POM_3))
    @test all(isapprox.(preferred_losses, (0.5, 0.5, 0.5)))

    zero_state = (N=0.0, POM_1=0.0, POM_2=0.0, POM_3=0.0, B_1=1.0)
    zero_args = (0.0, 0.0, 0.0, 0.0, Tuple(getproperty(zero_state, name) for name in names)...)
    @test all(iszero(bgc(Val(name), zero_args...)) for name in (:POM_1, :POM_2, :POM_3, :B_1, :N))
end

@testset "Preferential grazing shares consumer capacity across prey" begin
    function grazing_model(formulation; half_saturation=1.0)
        overrides = merge(
            food_web_parameter_overrides(),
            (
                maximum_predation_rate=[0.0, 1.0],
                holling_half_saturation=fill(half_saturation, 2),
                living_palatability_matrix=[0.0 0.0; 0.8 0.8],
                living_assimilation_matrix=ones(2, 2),
            ),
        )
        return construct(food_web_definition(; grazing=formulation); parameter_overrides=overrides)
    end
    prey_losses = (model, p, b) -> begin
        args = food_web_args(model, (P_1=p, B_1=b, Z_1=1.0))
        return (-model(Val(:P_1), args...), -model(Val(:B_1), args...))
    end

    expected_total = 0.8 / (1.0 + 0.8)
    proportional = grazing_model(PreferentialGrazing())
    @test sum(prey_losses(proportional, 1.0, 0.0)) ≈ expected_total
    @test sum(prey_losses(proportional, 0.5, 0.5)) ≈ expected_total

    switching = grazing_model(PreferentialGrazing(; switching_exponent=2))
    switched = prey_losses(switching, 0.75, 0.25)
    @test sum(switched) ≈ expected_total
    @test switched[1] / switched[2] ≈ 9.0

    zero_proportional = grazing_model(PreferentialGrazing(); half_saturation=0.0)
    zero_switching = grazing_model(PreferentialGrazing(; switching_exponent=2); half_saturation=0.0)
    @test prey_losses(zero_proportional, 0.0, 0.0) == (0.0, 0.0)
    @test prey_losses(zero_switching, 0.0, 0.0) == (0.0, 0.0)
end
