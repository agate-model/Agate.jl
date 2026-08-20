using Test

using Agate.Configuration: Population, Pool, realize_components
using Agate.Factories: default_components, default_processes
using Agate.Processes:
    Smith,
    Monod,
    Growth,
    Light,
    Nutrients,
    NutrientResponse,
    FixedStoichiometry,
    Grazing,
    ModelDefinition,
    ParameterSlot,
    ParameterRequirementIdentity,
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
    normalize_model,
    participants,
    process_id,
    process_kind,
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
    @test participants(grazing) == (
        consumer=(:Z,), resource=(:P, :B), unassimilated_destination=(:D,)
    )
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
    @test_throws ArgumentError Growth(; population=:P, factors=NamedTuple())
    @test_throws ArgumentError Growth(;
        population=:P,
        populations=(:P,),
        factors=(light=Light(:smith; driver=:PAR),),
    )

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
end

@testset "NiPiZD normalized process contract" begin
    factory = Agate.Models.NiPiZD.NiPiZDFactory()
    definition = ModelDefinition(factory)
    normalized = normalize_model(definition)

    @test normalized.components == default_components(factory)
    @test length(normalized.parameters) == 15
    @test driver_identities(normalized) == (:PAR,)
    @test keys(normalized.processes) == (
        :grazing_Z_on_P,
        :growth_P,
        :linear_mortality_P,
        :linear_mortality_Z,
        :quadratic_mortality_Z,
        :remineralization_D,
    )
    @test Tuple(process_id(process) for process in values(normalized.processes)) ==
        keys(normalized.processes)
    @test length(parameter_requirements(normalized)) == 18
    @test length(parameter_bindings(normalized)) == 18

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
        consumer=(:Z,), resource=(:P,), unassimilated_destination=(:D,)
    )
    @test participants(normalized.processes.remineralization_D) ==
        (source=(:D,), destination=(:N,))

    @test process_kind(normalized.processes.growth_P) === :growth
    @test formulation_tag(formulation(normalized.processes.grazing_Z_on_P.process)) ===
        :preferential
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
    ).axis_tracers == ((:P_1, :P_2),)
    @test application(:maximum_predation_rate, :grazing_Z_on_P, :maximum_rate).axis_tracers ==
        ((:Z_1, :Z_2),)
    @test application(
        :protection, :grazing_Z_on_P, :protection; path=(:palatability, :default)
    ).axis_tracers == ((:P_1, :P_2),)
    @test application(:palatability_matrix, :grazing_Z_on_P, :palatability).axis_tracers ==
        ((:Z_1, :Z_2), (:P_1, :P_2))
    remineralization_rate = application(
        :detritus_remineralization, :remineralization_D, :rate
    )
    @test remineralization_rate.binding.requirement.shape === :scalar
    @test remineralization_rate.axis_tracers == ((:D,),)
    @test isempty(
        application(
            :mortality_export_fraction,
            :linear_mortality_P,
            :export_fraction;
            path=(:routing,),
        ).axis_tracers,
    )

    linear_P_rate = only(
        filter(parameter_requirements(normalized.processes.linear_mortality_P)) do requirement
            requirement.identity.slot === :rate
        end,
    )
    @test parameter_name(normalized, linear_P_rate) === :linear_mortality

    processes = default_processes(factory)
    reversed_names = reverse(keys(processes))
    reversed_processes = NamedTuple{reversed_names}(reverse(values(processes)))
    reordered = normalize_model(
        ModelDefinition(;
            components=default_components(factory), processes=reversed_processes
        ),
    )
    @test keys(reordered.processes) == keys(normalized.processes)
    @test all(
        getproperty(reordered.processes, name).process ===
        getproperty(normalized.processes, name).process for name in keys(normalized.processes)
    )

    invalid = ModelDefinition(;
        components=default_components(factory),
        processes=(bad=Grazing(:preferential; consumer=:Z, resource=:missing),),
    )
    @test_throws ArgumentError normalize_model(invalid)

    @test_throws ArgumentError normalize_model(
        ModelDefinition(;
            components=default_components(factory),
            processes=default_processes(factory),
            parameters=(),
        ),
    )
    @test_throws ArgumentError normalize_model(
        ModelDefinition(;
            components=default_components(factory),
            processes=default_processes(factory),
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
