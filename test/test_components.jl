using Test
using Agate.Configuration:
    Population, Pool, currency, states, reference_state, variable_states, state_element,
    size_structure, realize_model_layout, component_classes, component_state_tracers,
    component_tracers, state_tracers, component_diameters
using Agate.ModelFamilies: default_components
@testset "Component authoring" begin
    population = Population(; states=:nitrogen, reference_state=:nitrogen,
        size_structure=(n=3, min_esd=1.0, max_esd=100.0, splitting=:log_splitting),
    )
    pool = Pool(:carbon)

    @test states(population) == (:nitrogen,)
    @test reference_state(population) === :nitrogen
    @test isempty(variable_states(population))
    @test state_element(population, :nitrogen) === :nitrogen
    @test size_structure(population).n == 3
    @test currency(pool) === :carbon
    @test isnothing(size_structure(pool))

    layout = realize_model_layout((P=population, D=pool); scalar_type=Float32)
    @test component_tracers(layout, :P) == (:P_1, :P_2, :P_3)
    @test collect(component_diameters(layout, :P)) ≈ Float32[1, 10, 100]
    @test isnothing(component_diameters(layout, :D))

    scalar_population = realize_model_layout((
        B=Population(; states=:carbon, reference_state=:carbon),
    ))
    @test component_tracers(scalar_population, :B) == (:B,)
    @test isnothing(component_diameters(scalar_population, :B))

    @test_throws ArgumentError Population(; states=nothing, reference_state=:carbon)
    @test_throws ArgumentError Pool(nothing)
    @test_throws ArgumentError realize_model_layout((P=population, P_1=Pool(:nitrogen)))
end

@testset "Multi-state population identity and realization" begin
    population = Population(;
        states=(:carbon, :nitrogen, :phosphorus),
        reference_state=:carbon,
        size_structure=[2.0, 10.0],
    )
    @test states(population) == (:carbon, :nitrogen, :phosphorus)
    @test reference_state(population) === :carbon
    @test variable_states(population) == (:nitrogen, :phosphorus)
    @test state_element(population, :nitrogen) === :nitrogen
    @test_throws ArgumentError state_element(population, :iron)

    layout = realize_model_layout((P=population, DIN=Pool(:nitrogen)); scalar_type=Float32)
    @test component_classes(layout, :P) == (:P_1, :P_2)
    @test length(component_classes(layout, :P)) == 2
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
    @test component_diameters(layout, :P) == (2.0f0, 10.0f0)
    grouped = realize_model_layout(
        (P=population, DIN=Pool(:nitrogen));
        scalar_type=Float32,
        population_groups=(P=(:P,),),
        group_diameters=(P=[2.0, 10.0],),
    )
    @test component_classes(grouped, :P) == component_classes(layout, :P)
    @test component_tracers(grouped, :P) == component_tracers(layout, :P)
    @test component_state_tracers(grouped, :P) == component_state_tracers(layout, :P)
    @test component_diameters(grouped, :P) == component_diameters(layout, :P)

    @test_throws ArgumentError Population(; states=(:carbon,), reference_state=nothing)
    @test_throws ArgumentError Population(; states=(), reference_state=:carbon)
    @test_throws ArgumentError Population(; states=(:carbon, :carbon), reference_state=:carbon)
    @test_throws ArgumentError Population(; states=(:carbon, 1), reference_state=:carbon)
    @test_throws ArgumentError Population(; states=(:nitrogen,), reference_state=:carbon)

    photoacclimating = Population(;
        states=(:carbon, :nitrogen, :chlorophyll),
        reference_state=:carbon,
    )
    @test variable_states(photoacclimating) == (:nitrogen, :chlorophyll)
    @test Tuple(state_element(photoacclimating, state) for state in states(photoacclimating)) ==
        (:carbon, :nitrogen, nothing)

    abundance = Population(; states=:abundance, reference_state=:abundance)
    @test isnothing(state_element(abundance, reference_state(abundance)))
end

@testset "Structured pool class layout" begin
    family = Agate.Models.NiPiZD.NiPiZDFamily()
    base_components = default_components(family)
    components = merge(
        base_components,
        (POM=Pool(:nitrogen; size_structure=[0.5, 5.0, 50.0]),),
    )
    layout = realize_model_layout(
        components;
        scalar_type=Float32,
        interaction_roles=(consumers=(:Z,), prey=(:P,)),
    )

    @test component_tracers(layout, :POM) == (:POM_1, :POM_2, :POM_3)
    @test component_diameters(layout, :POM) == (0.5f0, 5.0f0, 50.0f0)
end

@testset "NiPiZD component declaration" begin
    family = Agate.Models.NiPiZD.NiPiZDFamily()
    components = default_components(family)
    @test keys(components) == (:P, :Z, :N, :D)
    @test components.P isa Population
    @test components.Z isa Population
    @test all(
        state_element(getproperty(components, name), reference_state(getproperty(components, name))) === :nitrogen
        for name in (:P, :Z)
    )
    @test all(currency(getproperty(components, name)) === :nitrogen for name in (:N, :D))
end
