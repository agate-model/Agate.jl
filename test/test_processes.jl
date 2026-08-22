using Test
using ForwardDiff

using Agate.Configuration: Population, Pool, realize_components
using Agate.ModelFamilies: default_components, default_processes
using Agate.Parameters: NoDefault, ParameterDefinition, ParameterProvision
using Agate.Processes:
    Smith,
    Geider,
    Monod,
    Liebig,
    Frank,
    Q10,
    MultiplicativeFactors,
    Growth,
    Light,
    Nutrients,
    NutrientResponse,
    FixedStoichiometry,
    Grazing,
    Consumption,
    Mortality,
    DirectRouting,
    ModelDefinition,
    ParameterSlot,
    ParameterRequirementIdentity,
    ProductRouting,
    parameter_bindings,
    parameter_requirements,
    parameter_slots,
    parameter_slot_bindings,
    parameter_name,
    resolve_parameter_applicability,
    drivers,
    driver_identities,
    formulation,
    formulation_tag,
    factors,
    factor_value,
    normalize_model,
    participants,
    process_id,
    process_kind,
    process_rate,
    rate_axes

@testset "Process authoring and normalization" begin
    symbolic_light = Light(:smith; driver=:PAR)
    concrete_light = Light(Smith(); driver=:PAR)
    symbolic_response = NutrientResponse(:monod; resource=:N)
    concrete_response = NutrientResponse(Monod(); resource=:N)
    symbolic_growth = Growth(;
        population=:P,
        factors=(light=symbolic_light, nutrients=symbolic_response),
    )
    concrete_growth = Growth(;
        population=:P,
        factors=(nutrients=concrete_response, light=concrete_light),
    )

    @test typeof(formulation(symbolic_light)) === Smith
    @test typeof(formulation(symbolic_response)) === Monod
    frank_nutrients = Nutrients(
        :frank;
        responses=(nitrogen=NutrientResponse(:monod; resource=:N),),
        sharpness=25,
    )
    @test formulation_tag(formulation(frank_nutrients)) === :frank
    @test formulation(frank_nutrients).sharpness == 25
    @test factors(symbolic_growth) == factors(concrete_growth)
    @test keys(factors(symbolic_growth)) == (:light, :nutrients)
    @test participants(symbolic_growth) == participants(concrete_growth) == (population=(:P,),)
    @test drivers(symbolic_growth) == drivers(concrete_growth) == (light=:PAR,)
    @test rate_axes(symbolic_growth) == (:population,)

    grazing = Grazing(
        :preferential;
        consumers=:Z,
        resources=(:P, :B),
        unassimilated_destination=:D,
    )
    @test grazing isa Consumption
    @test process_kind(grazing) === :consumption
    @test participants(grazing) == (consumer=(:Z,), resource=(:P, :B))
    @test formulation(grazing.routing) isa DirectRouting
    @test grazing.routing.retained === :D
    @test rate_axes(grazing) == (:consumer, :resource)
    @test isempty(factors(grazing))

    shared_driver_model = normalize_model(
        ModelDefinition(;
            components=(
                P=Population(; currency=:nitrogen), Z=Population(; currency=:nitrogen)
            ),
            processes=(
                growth_Z=Growth(;
                    population=:Z, factors=(light=Light(:smith; driver=:PAR),)
                ),
                growth_P=Growth(;
                    population=:P, factors=(light=Light(:smith; driver=:PAR),)
                ),
            ),
        ),
    )
    @test driver_identities(shared_driver_model) == (:PAR,)

    invalid_growth = ModelDefinition(;
        components=(P=Population(; currency=:nitrogen),),
        processes=(growth=Growth(;
            population=:P,
            factors=(nutrients=NutrientResponse(:monod; resource=:missing),),
        ),),
    )
    @test_throws ArgumentError normalize_model(invalid_growth)

    @test_throws ArgumentError Light(:unknown; driver=:PAR)
    @test_throws ArgumentError Nutrients(
        :liebig;
        responses=(nitrogen=NutrientResponse(:monod; resource=:N),),
        sharpness=25,
    )
    @test_throws ArgumentError Growth(; population=:P, factors=NamedTuple())
    @test_throws ArgumentError Growth(;
        population=:P,
        populations=(:P,),
        factors=(light=Light(:smith; driver=:PAR),),
    )
    @test_throws ArgumentError Grazing(
        :preferential;
        consumer=:Z,
        resource=:P,
        unassimilated_destination=:D,
        routing=ProductRouting(:direct; destination=:N),
    )

    redundant_growth_source = ModelDefinition(;
        components=(
            P=Population(; currency=:nitrogen),
            N=Pool(:nitrogen),
            D=Pool(:nitrogen),
        ),
        processes=(growth=Growth(;
            population=:P,
            source=:D,
            factors=(nutrients=NutrientResponse(:monod; resource=:N),),
        ),),
    )
    @test_throws ArgumentError normalize_model(redundant_growth_source)

    @test parameter_slots(Smith()) == (
        ParameterSlot(:maximum_rate, (:population,)),
        ParameterSlot(:alpha, (:population,)),
    )
    @test parameter_slots(Monod()) == (
        ParameterSlot(:K, (:population,); qualify=:resource),
    )

    wrong_currency = ModelDefinition(;
        components=(
            P=Population(; currency=:carbon),
            DIC=Pool(:carbon),
            DIN=Pool(:phosphorus),
        ),
        processes=(growth=Growth(;
            population=:P,
            source=:DIC,
            factors=(
                light=Light(:geider; driver=:PAR),
                nutrients=Nutrients(
                    :liebig;
                    responses=(nitrogen=NutrientResponse(:monod; resource=:DIN),),
                ),
            ),
            stoichiometry=FixedStoichiometry(; reference=:carbon),
        ),),
    )
    @test_throws ArgumentError normalize_model(wrong_currency)

    # Normalization remains the trusted boundary when field constructors bypass authoring checks.
    bypassed_factor = ModelDefinition(;
        components=(P=Population(; currency=:nitrogen), N=Pool(:nitrogen)),
        processes=(growth=Growth(;
            population=:P,
            factors=(light=Light(Monod(), :PAR), nutrients=NutrientResponse(:monod; resource=:N)),
        ),),
    )
    @test_throws ArgumentError normalize_model(bypassed_factor)

    bypassed_mortality = ModelDefinition(;
        components=(P=Population(; currency=:nitrogen),),
        processes=(mortality=Mortality(Monod(), (:P,), nothing),),
    )
    @test_throws ArgumentError normalize_model(bypassed_mortality)
end

@testset "Compositional factor mathematics" begin
    biomass = 0.4
    maximum_rate = 2e-5
    light = 100.0
    alpha = 2e-6

    base = process_rate(MultiplicativeFactors(), biomass, maximum_rate)
    smith = factor_value(Smith(), light, maximum_rate, alpha)
    monod = factor_value(Monod(), 0.5, 0.2)
    @test isapprox(
        base * smith * monod,
        maximum_rate * biomass *
        (alpha * light / sqrt(maximum_rate^2 + (alpha * light)^2)) *
        (0.5 / 0.7);
        rtol=1e-14,
    )

    chlorophyll_to_carbon = 0.02
    geider = factor_value(
        Geider(), light, maximum_rate, alpha, chlorophyll_to_carbon
    )
    @test geider == Agate.Library.Photosynthesis.geider_light_response(
        light, alpha, maximum_rate, chlorophyll_to_carbon
    )
    @test factor_value(Geider(), light, 0.0, alpha, chlorophyll_to_carbon) == 0.0
    limitations = (factor_value(Monod(), 3.0, 0.5), factor_value(Monod(), 0.2, 0.5))
    liebig = factor_value(Liebig(), limitations)
    frank = factor_value(Frank(25), limitations)
    @test frank ≈ Agate.Library.Nutrients.FrankTNorm(25)(limitations)
    frank_gradient = ForwardDiff.gradient(
        x -> factor_value(Frank(25), (x[1], x[2])), collect(limitations)
    )
    @test all(isfinite, frank_gradient)
    temperature = factor_value(Q10(), 30.0, 2.0, 20.0)
    @test temperature == Agate.Library.Temperature.q10_temperature_factor(30.0, 2.0, 20.0)
    expected_liebig = min(3.0 / 3.5, 0.2 / 0.7)
    expected_geider = 1 - exp(-alpha * chlorophyll_to_carbon * light / maximum_rate)
    @test isapprox(
        base * geider * liebig * temperature,
        maximum_rate * biomass * expected_geider * expected_liebig * 2;
        rtol=1e-14,
    )
end

@testset "NiPiZD normalized process contract" begin
    family = Agate.Models.NiPiZD.NiPiZDFamily()
    definition = ModelDefinition(family)
    normalized = normalize_model(definition)

    @test normalized.components == default_components(family)
    @test length(normalized.parameters) == 15
    @test driver_identities(normalized) == (:PAR,)
    @test keys(normalized.processes) == (
        :grazing_Z_on_P,
        :growth_P,
        :linear_mortality_P_to_D,
        :linear_mortality_P_to_N,
        :linear_mortality_Z_to_N,
        :quadratic_mortality_Z_to_D,
        :remineralization_D,
    )
    @test Tuple(process_id(process) for process in values(normalized.processes)) ==
        keys(normalized.processes)
    @test length(parameter_requirements(normalized)) == 16
    @test length(parameter_bindings(normalized)) == 16

    growth = normalized.processes.growth_P
    light_slots = parameter_slot_bindings(
        normalized, growth, (:factors, :light), growth.process.factors.light.formulation
    )
    @test keys(light_slots) == Tuple(slot.name for slot in parameter_slots(Smith()))
    @test Tuple(binding.requirement.identity.slot for binding in values(light_slots)) ==
        keys(light_slots)
    @test Tuple(binding.parameter for binding in values(light_slots)) ==
        (:maximum_growth_rate, :alpha)

    grazing = normalized.processes.grazing_Z_on_P
    grazing_slots = parameter_slot_bindings(
        normalized, grazing, (), grazing.process.formulation
    )
    @test (
        grazing_slots.palatability.runtime_path, grazing_slots.assimilation.runtime_path
    ) == ((:interactions, :palatability), (:interactions, :assimilation))

    @test participants(normalized.processes.growth_P) == (population=(:P,),)
    @test drivers(normalized.processes.growth_P) == (light=:PAR,)
    @test participants(normalized.processes.grazing_Z_on_P) == (
        consumer=(:Z,), resource=(:P,)
    )
    @test process_kind(normalized.processes.grazing_Z_on_P) === :consumption
    @test participants(normalized.processes.remineralization_D) ==
        (source=(:D,), destination=(:N,))

    @test process_kind(normalized.processes.growth_P) === :growth
    @test rate_axes(normalized.processes.growth_P) == (:population,)
    @test rate_axes(normalized.processes.grazing_Z_on_P) == (:consumer, :resource)

    layout = realize_components(normalized.components)
    applicability = resolve_parameter_applicability(normalized, layout)
    function application(parameter, process, slot; path=())
        return only(
            filter(applicability) do item
                identity = item.binding.requirement.identity
                item.binding.parameter === parameter && identity.process === process &&
                    identity.path == path && identity.slot === slot
            end,
        )
    end
    @test application(
        :maximum_growth_rate, :growth_P, :maximum_rate; path=(:factors, :light)
    ).axis_classes == ((:P_1, :P_2),)
    @test application(:maximum_predation_rate, :grazing_Z_on_P, :maximum_rate).axis_classes ==
        ((:Z_1, :Z_2),)
    @test application(
        :protection, :grazing_Z_on_P, :protection; path=(:palatability, :default)
    ).axis_classes == ((:P_1, :P_2),)
    @test application(:palatability_matrix, :grazing_Z_on_P, :palatability).axis_classes ==
        ((:Z_1, :Z_2), (:P_1, :P_2))
    remineralization_rate = application(
        :detritus_remineralization, :remineralization_D, :rate
    )
    @test remineralization_rate.binding.requirement.shape === :scalar
    @test remineralization_rate.axis_classes == ((:D,),)
    linear_P_to_N_rate = only(
        filter(parameter_requirements(normalized.processes.linear_mortality_P_to_N)) do requirement
            requirement.identity.slot === :rate
        end,
    )
    linear_P_to_D_rate = only(
        filter(parameter_requirements(normalized.processes.linear_mortality_P_to_D)) do requirement
            requirement.identity.slot === :rate
        end,
    )
    @test parameter_name(normalized, linear_P_to_N_rate) === :linear_mortality
    @test parameter_name(normalized, linear_P_to_D_rate) === :linear_detrital_mortality

    processes = default_processes(family)
    reversed_names = reverse(keys(processes))
    reversed_processes = NamedTuple{reversed_names}(reverse(values(processes)))
    reordered = normalize_model(
        ModelDefinition(;
            components=default_components(family), processes=reversed_processes
        ),
    )
    @test keys(reordered.processes) == keys(normalized.processes)
    @test all(
        getproperty(reordered.processes, name).process ===
        getproperty(normalized.processes, name).process for name in keys(normalized.processes)
    )

    invalid = ModelDefinition(;
        components=default_components(family),
        processes=(bad=Grazing(:preferential; consumer=:Z, resource=:missing),),
    )
    @test_throws ArgumentError normalize_model(invalid)

    @test_throws ArgumentError normalize_model(
        ModelDefinition(;
            components=default_components(family),
            processes=default_processes(family),
            parameters=(),
        ),
    )
    @test_throws ArgumentError normalize_model(
        ModelDefinition(;
            components=default_components(family),
            processes=default_processes(family),
            parameters=(definition.parameters..., first(definition.parameters)),
        ),
    )
end

@testset "Parameter requirement identity" begin
    a = ParameterRequirementIdentity(
        :growth_P,
        (:factors, :nutrients),
        Monod(),
        :K;
        qualifier=(resource=:N, population=:P),
    )
    b = ParameterRequirementIdentity(
        :growth_P,
        (:factors, :nutrients),
        :monod,
        :K;
        qualifier=(population=:P, resource=:N),
    )

    @test a == b
    @test a != ParameterRequirementIdentity(
        :other_growth, (:factors, :nutrients), :monod, :K; qualifier=(resource=:N, population=:P)
    )
end

@testset "Parameter provision inference" begin
    components = (P=Population(; currency=:nitrogen, size_structure=[1.0]),)
    processes = (
        growth=Growth(;
            population=:P,
            factors=(
                light_a=Light(:smith; driver=:PAR),
                light_b=Light(:smith; driver=:PAR),
            ),
        ),
    )
    parameter(name, slot, path) = ParameterDefinition(
        name, NoDefault();
        shape=:vector,
        provides=ParameterProvision(:growth, slot; path),
    )
    parameters(max_a) = (
        ParameterDefinition(:max_a, NoDefault(); shape=:vector, provides=max_a),
        parameter(:alpha_a, :alpha, (:factors, :light_a)),
        parameter(:max_b, :maximum_rate, (:factors, :light_b)),
        parameter(:alpha_b, :alpha, (:factors, :light_b)),
    )
    definition(parameters) = ModelDefinition(; components, processes, parameters)

    normalized = normalize_model(definition(parameters(
        ParameterProvision(:growth, :maximum_rate; path=(:factors, :light_a))
    )))
    max_a_requirement = only(filter(parameter_requirements(normalized)) do requirement
        identity = requirement.identity
        identity.path == (:factors, :light_a) && identity.slot === :maximum_rate
    end)
    @test max_a_requirement.identity.formulation === :smith
    @test parameter_name(normalized, max_a_requirement) === :max_a

    @test_throws ArgumentError normalize_model(
        definition(parameters(ParameterProvision(:growth, :maximum_rate)))
    )
    @test_throws ArgumentError normalize_model(
        definition(parameters(ParameterProvision(:growth, :missing)))
    )
end
