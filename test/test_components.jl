using Test
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers

using Agate.Configuration:
    Population, Pool, currency, states, state_currency, size_structure, realize_components,
    realize_component_groups, component_classes, component_state_tracers, component_tracers,
    component_indices, state_tracers, state_indices,
    component_diameters, component_class_count, parse_community
using Agate.ModelFamilies: default_components
using Agate.Parameters: ParameterProvision, ParameterDefinition, ConstDefault, NoDefault
using Agate.Processes:
    ModelDefinition, Mortality, Remineralization, normalize_model, resolve_parameter_applicability

@testset "Component authoring" begin
    population = Population(;
        currency=:nitrogen,
        size_structure=(n=3, min_esd=1.0, max_esd=100.0, splitting=:log_splitting),
    )
    pool = Pool(:carbon)

    @test currency(population) === :nitrogen
    @test states(population) == (nitrogen=:nitrogen,)
    @test state_currency(population, :nitrogen) === :nitrogen
    @test size_structure(population).n == 3
    @test currency(pool) === :carbon
    @test isnothing(size_structure(pool))

    layout = realize_components((P=population, D=pool); scalar_type=Float32)
    @test layout.tracer_order == (:P_1, :P_2, :P_3, :D)
    @test layout.component_classes == (P=(:P_1, :P_2, :P_3), D=(:D,))
    @test layout.component_state_tracers == (P=(nitrogen=(:P_1, :P_2, :P_3),), D=nothing)
    @test layout.component_tracers == (P=(:P_1, :P_2, :P_3), D=(:D,))
    @test layout.component_indices == (P=(1, 2, 3), D=(4,))
    @test collect(layout.component_diameters.P) ≈ Float32[1, 10, 100]
    @test isnothing(layout.component_diameters.D)

    scalar_population = realize_components((B=Population(; currency=:carbon),))
    @test scalar_population.tracer_order == (:B,)
    @test isnothing(scalar_population.component_diameters.B)

    structured_pool = Pool(:carbon; size_structure=[1.0, 10.0, 100.0])
    generic = realize_components((P=population, POM=structured_pool); scalar_type=Float32)
    @test generic.tracer_order == (:P_1, :P_2, :P_3, :POM_1, :POM_2, :POM_3)
    @test component_tracers(generic, :POM) == (:POM_1, :POM_2, :POM_3)
    @test component_indices(generic, :POM) == (4, 5, 6)
    @test component_diameters(generic, :POM) == (1.0f0, 10.0f0, 100.0f0)
    @test component_class_count(generic, :POM) == 3
    @test_throws ArgumentError Population(; currency=nothing)
    @test_throws ArgumentError Pool(nothing)
    @test_throws ArgumentError realize_components((P=population, P_1=Pool(:nitrogen)))
end

@testset "Multi-state population identity and realization" begin
    population = Population(;
        states=(carbon=:carbon, nitrogen=:nitrogen, phosphorus=:phosphorus),
        size_structure=[2.0, 10.0],
    )
    @test states(population) == (carbon=:carbon, nitrogen=:nitrogen, phosphorus=:phosphorus)
    @test state_currency(population, :nitrogen) === :nitrogen
    @test_throws ArgumentError currency(population)
    @test_throws ArgumentError state_currency(population, :iron)

    layout = realize_components((P=population, DIN=Pool(:nitrogen)); scalar_type=Float32)
    @test component_classes(layout, :P) == (:P_1, :P_2)
    @test component_class_count(layout, :P) == 2
    @test component_tracers(layout, :P) == (
        :P_1_carbon, :P_1_nitrogen, :P_1_phosphorus,
        :P_2_carbon, :P_2_nitrogen, :P_2_phosphorus,
    )
    @test component_state_tracers(layout, :P) == (
        carbon=(:P_1_carbon, :P_2_carbon),
        nitrogen=(:P_1_nitrogen, :P_2_nitrogen),
        phosphorus=(:P_1_phosphorus, :P_2_phosphorus),
    )
    @test state_tracers(layout, :P, :nitrogen) == (:P_1_nitrogen, :P_2_nitrogen)
    @test state_indices(layout, :P, :carbon) == (1, 4)
    @test state_indices(layout, :P, :phosphorus) == (3, 6)
    @test component_diameters(layout, :P) == (2.0f0, 10.0f0)
    @test layout.tracer_order == (
        :P_1_carbon, :P_1_nitrogen, :P_1_phosphorus,
        :P_2_carbon, :P_2_nitrogen, :P_2_phosphorus, :DIN,
    )

    community = (P=(diameters=[2.0, 10.0], pft=Agate.Configuration.PFTSpecification()),)
    context = parse_community(Float32, community; biogeochem_tracers=(:DIN,))
    grouped = realize_component_groups((P=population, DIN=Pool(:nitrogen)), (P=(:P,),), context)
    @test context.class_symbols == [:P_1, :P_2]
    @test context.n_total == 2
    @test context.consumer_indices == [1, 2]
    @test context.prey_indices == [1, 2]
    @test component_classes(grouped, :P) == (:P_1, :P_2)
    @test grouped.tracer_order == (
        :DIN,
        :P_1_carbon, :P_1_nitrogen, :P_1_phosphorus,
        :P_2_carbon, :P_2_nitrogen, :P_2_phosphorus,
    )

    mortality = Mortality(:linear; population=:P)
    parameter = ParameterDefinition(
        :mortality_rate,
        NoDefault();
        shape=:vector,
        provides=ParameterProvision(:mortality_P, :rate),
    )
    definition = normalize_model(ModelDefinition(;
        components=(P=population,),
        processes=(mortality_P=mortality,),
        parameters=(parameter,),
    ))
    applicability = only(resolve_parameter_applicability(definition, realize_components((P=population,))))
    @test applicability.axis_components == ((:P,),)
    @test applicability.axis_classes == ((:P_1, :P_2),)

    @test_throws ArgumentError Population()
    @test_throws ArgumentError Population(; currency=:carbon, states=(carbon=:carbon,))
    @test_throws ArgumentError Population(; states=NamedTuple())
    @test_throws ArgumentError Population(; states=(carbon=nothing,))
    @test_throws ArgumentError Population(; states=(carbon=1,))
end

@testset "Structured pool class layout and parameter applicability" begin
    family = Agate.Models.NiPiZD.NiPiZDFamily()
    base_components = default_components(family)
    components = merge(
        base_components,
        (POM=Pool(:nitrogen; size_structure=[0.5, 5.0, 50.0]),),
    )
    community = default_nipizd_community()
    context = parse_community(
        Float32, community; biogeochem_tracers=(:N, :D, :POM_1, :POM_2, :POM_3)
    )
    layout = realize_component_groups(
        components, (P=(:P,), Z=(:Z,)), context
    )

    @test layout.tracer_order ==
          (:N, :D, :POM_1, :POM_2, :POM_3, :Z_1, :Z_2, :P_1, :P_2)
    @test component_indices(layout, :POM) == (3, 4, 5)
    @test component_diameters(layout, :POM) == (0.5f0, 5.0f0, 50.0f0)

    process = Remineralization(:linear; source=:POM, destination=:N)
    parameter = ParameterDefinition(
        :pom_remineralization,
        ConstDefault(0.1);
        provides=ParameterProvision(:remineralization_POM, :rate),
    )
    definition = normalize_model(
        ModelDefinition(;
            components,
            processes=(remineralization_POM=process,),
            parameters=(parameter,),
        ),
    )
    applicability = only(resolve_parameter_applicability(definition, layout))
    @test applicability.axis_components == ((:POM,),)
    @test applicability.axis_classes == ((:POM_1, :POM_2, :POM_3),)
end

@testset "NiPiZD component declaration" begin
    family = Agate.Models.NiPiZD.NiPiZDFamily()
    components = default_components(family)
    layout = realize_components(components)

    @test keys(components) == (:P, :Z, :N, :D)
    @test components.P isa Population
    @test components.Z isa Population
    @test components.N isa Pool
    @test components.D isa Pool
    @test all(currency(getproperty(components, name)) === :nitrogen for name in keys(components))

    @test layout.component_tracers == (
        P=(:P_1, :P_2), Z=(:Z_1, :Z_2), N=(:N,), D=(:D,)
    )
    @test collect(layout.component_diameters.P) ≈ [2.0, 10.0]
    @test collect(layout.component_diameters.Z) ≈ [20.0, 100.0]

    bgc = Agate.Models.NiPiZD.construct(; grid=dummy_grid(Float32))
    runtime_tracers = required_biogeochemical_tracers(bgc)
    @test Set(layout.tracer_order) == Set(runtime_tracers)
    @test runtime_tracers == (:N, :D, :Z_1, :Z_2, :P_1, :P_2)
end
