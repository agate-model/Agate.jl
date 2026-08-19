using JSON
using Test

using Agate.Configuration: Population, realize_components
using Agate.Factories: default_components, default_processes
using Agate.Processes:
    Smith,
    Monod,
    Growth,
    NutrientResponse,
    Grazing,
    ModelDefinition,
    ParameterRequirementIdentity,
    parameter_bindings,
    parameter_requirements,
    parameter_name,
    resolve_parameter_applicability,
    drivers,
    driver_identities,
    formulation,
    formulation_tag,
    normalize_model,
    participants,
    process_id,
    process_kind,
    rate_axes

@testset "Process authoring and normalization" begin
    symbolic_response = NutrientResponse(:monod; resource=:N)
    concrete_response = NutrientResponse(Monod(); resource=:N)
    symbolic_growth = Growth(
        :smith; population=:P, light=:PAR, limitation=symbolic_response
    )
    concrete_growth = Growth(
        Smith(); population=:P, light=:PAR, limitation=concrete_response
    )

    @test typeof(formulation(symbolic_growth)) === Smith
    @test typeof(formulation(symbolic_response)) === Monod
    @test formulation_tag(formulation(symbolic_growth)) ==
        formulation_tag(formulation(concrete_growth))
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

    shared_driver_model = normalize_model(
        ModelDefinition(;
            components=(
                P=Population(; currency=:nitrogen), Z=Population(; currency=:nitrogen)
            ),
            processes=(
                growth_Z=Growth(:smith; population=:Z, light=:PAR),
                growth_P=Growth(:smith; population=:P, light=:PAR),
            ),
        ),
    )
    @test driver_identities(shared_driver_model) == (:PAR,)

    @test_throws ArgumentError Growth(:unknown; population=:P, light=:PAR)
    @test_throws ArgumentError Growth(:smith; population=:P)
    @test_throws ArgumentError Growth(
        :smith; population=:P, populations=(:P,), light=:PAR
    )
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

    @test participants(normalized.processes.growth_P) == (population=(:P,),)
    @test drivers(normalized.processes.growth_P) == (light=:PAR,)
    @test participants(normalized.processes.grazing_Z_on_P) == (
        consumer=(:Z,), resource=(:P,), unassimilated_destination=(:D,)
    )
    @test participants(normalized.processes.remineralization_D) ==
        (source=(:D,), destination=(:N,))

    target = JSON.parsefile(
        joinpath(@__DIR__, "reference", "nipizd_model_recipe_v3_target.json")
    )["recipe"]
    target_processes = target["processes"]
    @test Set(String.(keys(normalized.processes))) == Set(keys(target_processes))
    for (name, process) in pairs(normalized.processes)
        expected = target_processes[String(name)]
        @test String(process_kind(process)) == expected["kind"]
        @test String(formulation_tag(formulation(process))) == expected["formulation"]
        @test String.(collect(rate_axes(process))) == expected["rate_axes"]
    end

    for (parameter, specification) in pairs(target["parameter_bindings"])
        for expected in specification["provides"]
            matches = filter(parameter_bindings(normalized)) do binding
                identity = binding.requirement.identity
                qualifier = Dict(String(k) => String(v) for (k, v) in pairs(identity.qualifier))
                binding.parameter === Symbol(parameter) &&
                    identity.process === Symbol(expected["process"]) &&
                    identity.path == Tuple(Symbol.(expected["path"])) &&
                    identity.formulation === Symbol(expected["formulation"]) &&
                    identity.slot === Symbol(expected["slot"]) &&
                    qualifier == expected["qualifier"] &&
                    String.(collect(binding.requirement.axes)) == expected["axes"]
            end
            @test length(matches) == 1
        end
    end

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
    @test application(:maximum_growth_rate, :growth_P, :maximum_rate).axis_tracers ==
        ((:P_1, :P_2),)
    @test application(:maximum_predation_rate, :grazing_Z_on_P, :maximum_rate).axis_tracers ==
        ((:Z_1, :Z_2),)
    @test application(
        :protection, :grazing_Z_on_P, :protection; path=(:palatability, :default)
    ).axis_tracers == ((:P_1, :P_2),)
    @test application(:palatability_matrix, :grazing_Z_on_P, :palatability).axis_tracers ==
        ((:Z_1, :Z_2), (:P_1, :P_2))
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
        (:limitation,),
        Monod(),
        :K;
        qualifier=(resource=:N, population=:P),
    )
    b = ParameterRequirementIdentity(
        :growth_P,
        (:limitation,),
        :monod,
        :K;
        qualifier=(population=:P, resource=:N),
    )

    @test a == b
    @test a != ParameterRequirementIdentity(
        :other_growth, (:limitation,), :monod, :K; qualifier=(resource=:N, population=:P)
    )
end
