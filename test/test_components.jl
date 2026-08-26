using Test
using Agate.Configuration:
    Population, Pool, currency, states, state_currency, size_structure, realize_components,
    realize_model_layout, component_classes, component_class_indices, component_state_tracers,
    component_tracers, component_indices, state_tracers, state_indices,
    component_diameters, component_class_count
using Agate.ModelFamilies: default_components
using Agate.Parameters: Parameter, ConstantDefault, NoDefault
using Agate.Processes:
    ModelDefinition, Mortality, Remineralization, LinearMortality, LinearRemineralization,
    normalize_model, resolve_parameter_applicability

@testset "Component authoring" begin
    population = Population(:nitrogen;
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
    @test layout.tracer_order == (:D, :P_1, :P_2, :P_3)
    @test component_tracers(layout, :P) == (:P_1, :P_2, :P_3)
    @test component_class_indices(layout, :P) == (1, 2, 3)
    @test collect(layout.component_diameters.P) ≈ Float32[1, 10, 100]
    @test isnothing(layout.component_diameters.D)

    scalar_population = realize_components((B=Population(:carbon),))
    @test scalar_population.tracer_order == (:B,)
    @test isnothing(scalar_population.component_diameters.B)

    @test_throws ArgumentError Population(nothing)
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
    @test component_class_indices(layout, :P) == (1, 2)
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
    for (state, expected) in pairs((carbon=(2, 5), nitrogen=(3, 6), phosphorus=(4, 7)))
        @test state_indices(layout, :P, state) == expected
    end
    @test component_diameters(layout, :P) == (2.0f0, 10.0f0)
    @test layout.tracer_order == (
        :DIN,
        :P_1_carbon, :P_1_nitrogen, :P_1_phosphorus,
        :P_2_carbon, :P_2_nitrogen, :P_2_phosphorus,
    )

    grouped = realize_model_layout(
        (P=population, DIN=Pool(:nitrogen));
        scalar_type=Float32,
        population_groups=(P=(:P,),),
        group_diameters=(P=[2.0, 10.0],),
    )
    @test grouped.tracer_order == layout.tracer_order
    @test grouped.component_classes == layout.component_classes
    @test grouped.component_state_indices == layout.component_state_indices

    mortality = Mortality(
        LinearMortality(); populations=:P, bindings=(rate=:mortality_rate,)
    )
    parameter = Parameter(NoDefault())
    definition = normalize_model(ModelDefinition(;
        components=(P=population,),
        processes=(mortality_P=mortality,),
        parameters=(mortality_rate=parameter,),
    ))
    applicability = only(resolve_parameter_applicability(definition, realize_components((P=population,))))
    @test applicability.axis_components == ((:P,),)
    @test applicability.axis_classes == ((:P_1, :P_2),)

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
    realization = default_nipizd_realization()
    layout = realize_model_layout(
        components;
        scalar_type=Float32,
        realization...,
        interaction_roles=(consumers=(:Z,), prey=(:P,)),
    )

    @test layout.tracer_order ==
          (:N, :D, :POM_1, :POM_2, :POM_3, :Z_1, :Z_2, :P_1, :P_2)
    @test component_indices(layout, :POM) == (3, 4, 5)
    @test component_diameters(layout, :POM) == (0.5f0, 5.0f0, 50.0f0)

    process = Remineralization(
        LinearRemineralization();
        sources=:POM,
        destination=:N,
        bindings=(rate=:pom_remineralization,),
    )
    parameter = Parameter(ConstantDefault(0.1))
    definition = normalize_model(
        ModelDefinition(;
            components,
            processes=(remineralization_POM=process,),
            parameters=(pom_remineralization=parameter,),
        ),
    )
    applicability = only(resolve_parameter_applicability(definition, layout))
    @test applicability.axis_components == ((:POM,),)
    @test applicability.axis_classes == ((:POM_1, :POM_2, :POM_3),)
end

@testset "Direct and family population realization share layout facts" begin
    components = (
        N=Pool(:nitrogen), D=Pool(:nitrogen),
        Z=Population(:nitrogen; size_structure=[20.0, 100.0]),
        P=Population(:nitrogen; size_structure=[1.0, 10.0]),
    )
    direct = realize_model_layout(
        components; interaction_roles=(consumers=(:Z,), prey=(:P,))
    )
    family = realize_model_layout(
        components;
        population_groups=(Z=(:Z,), P=(:P,)),
        group_diameters=(Z=[20.0, 100.0], P=[1.0, 10.0]),
        interaction_roles=(consumers=(:Z,), prey=(:P,)),
    )
    for field in (
        :tracer_order, :component_classes, :component_state_indices, :component_indices,
        :class_symbols, :group_indices, :diameters, :consumer_indices, :prey_indices,
    )
        @test getfield(direct, field) == getfield(family, field)
    end
end

@testset "Structured subgroup indices" begin
    layout = realize_model_layout(
        (P=Population(:carbon), Z=Population(:carbon));
        population_groups=(P=(:diat, :dino), Z=(:micro, :meso)),
        group_diameters=(micro=[20.0, 30.0], meso=[100.0], diat=[2.0, 5.0], dino=[10.0]),
        interaction_roles=(consumers=(:micro, :meso), prey=(:diat, :dino)),
    )
    @test layout.group_indices == (
        micro=(1, 2), meso=(3,), diat=(4, 5), dino=(6,),
    )
    @test layout.group_local_indices == (1, 2, 1, 1, 2, 1)
    @test layout.component_local_indices == (1, 2, 3, 1, 2, 3)
    @test (layout.consumer_indices, layout.prey_indices) == ((1, 2, 3), (4, 5, 6))
end

@testset "NiPiZD component declaration" begin
    family = Agate.Models.NiPiZD.NiPiZDFamily()
    components = default_components(family)
    @test keys(components) == (:P, :Z, :N, :D)
    @test components.P isa Population
    @test components.Z isa Population
    @test all(currency(getproperty(components, name)) === :nitrogen for name in keys(components))
end
