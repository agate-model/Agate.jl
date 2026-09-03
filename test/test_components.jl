using Test
using Agate.Components:
    Plankton, Pool, element, states, reference_state, variable_states, state_element,
    size_structure, realize_model_layout, component_entities, component_state_tracers,
    component_tracers, state_tracers, component_diameters
using Agate.ModelFamilies: default_components
@testset "Component authoring" begin
    plankton = Plankton(; states=(nitrogen=:nitrogen,), reference_state=:nitrogen,
        size_structure=(n=3, min_esd=1.0, max_esd=100.0, spacing=:log),
    )
    pool = Pool(:carbon)

    @test states(plankton) == (:nitrogen,)
    @test reference_state(plankton) === :nitrogen
    @test isempty(variable_states(plankton))
    @test state_element(plankton, :nitrogen) === :nitrogen
    @test size_structure(plankton).n == 3
    @test element(pool) === :carbon

    layout = realize_model_layout((P=plankton, D=pool); scalar_type=Float32)
    @test component_tracers(layout, :P) == (:P_1, :P_2, :P_3)
    @test component_entities(layout, :D) == (:D,)
    @test component_tracers(layout, :D) == (:D,)
    @test collect(component_diameters(layout, :P)) ≈ Float32[1, 10, 100]
    @test isnothing(component_diameters(layout, :D))

    scalar_plankton_layout = realize_model_layout((
        B=Plankton(; states=(carbon=:carbon,), reference_state=:carbon),
    ))
    @test component_tracers(scalar_plankton_layout, :B) == (:B,)
    @test isnothing(component_diameters(scalar_plankton_layout, :B))

    @test_throws ArgumentError Plankton(; states=nothing, reference_state=:carbon)
    @test_throws MethodError Pool(nothing)
    @test_throws ArgumentError realize_model_layout((P=plankton, P_1=Pool(:nitrogen)))
end

@testset "Multi-state plankton identity and realization" begin
    plankton = Plankton(;
        states=(biomass=:carbon, reserve_n=:nitrogen, reserve_p=:phosphorus),
        reference_state=:biomass,
        size_structure=[2.0, 10.0],
    )
    @test states(plankton) == (:biomass, :reserve_n, :reserve_p)
    @test reference_state(plankton) === :biomass
    @test variable_states(plankton) == (:reserve_n, :reserve_p)
    @test state_element(plankton, :reserve_n) === :nitrogen
    @test_throws ArgumentError state_element(plankton, :iron)

    layout = realize_model_layout((P=plankton, DIN=Pool(:nitrogen)); scalar_type=Float32)
    @test component_entities(layout, :P) == (:P_1, :P_2)
    @test component_tracers(layout, :P) == (
        :P_1_biomass, :P_1_reserve_n, :P_1_reserve_p,
        :P_2_biomass, :P_2_reserve_n, :P_2_reserve_p,
    )
    @test component_state_tracers(layout, :P) == (
        biomass=(:P_1_biomass, :P_2_biomass),
        reserve_n=(:P_1_reserve_n, :P_2_reserve_n),
        reserve_p=(:P_1_reserve_p, :P_2_reserve_p),
    )
    @test state_tracers(layout, :P, :reserve_n) == (:P_1_reserve_n, :P_2_reserve_n)
    @test component_diameters(layout, :P) == (2.0f0, 10.0f0)
    explicit_pft_layout = realize_model_layout(
        (P=plankton, DIN=Pool(:nitrogen));
        scalar_type=Float32,
        plankton_pfts=(P=(P=[2.0, 10.0],),),
    )
    @test component_entities(explicit_pft_layout, :P) == component_entities(layout, :P)
    @test component_tracers(explicit_pft_layout, :P) == component_tracers(layout, :P)
    @test component_state_tracers(explicit_pft_layout, :P) == component_state_tracers(layout, :P)
    @test component_diameters(explicit_pft_layout, :P) == component_diameters(layout, :P)

    @test_throws ArgumentError Plankton(; states=(carbon=:carbon,), reference_state=nothing)
    @test_throws ArgumentError Plankton(; states=NamedTuple(), reference_state=:carbon)
    @test_throws ArgumentError Plankton(; states=(biomass=1,), reference_state=:biomass)
    @test_throws ArgumentError Plankton(; states=(a=:carbon, b=:carbon), reference_state=:a)
    @test_throws ArgumentError Plankton(; states=(nitrogen=:nitrogen,), reference_state=:carbon)

    photoacclimating = Plankton(;
        states=(carbon=:carbon, nitrogen=:nitrogen, chlorophyll=nothing),
        reference_state=:carbon,
    )
    @test variable_states(photoacclimating) == (:nitrogen, :chlorophyll)
    @test Tuple(state_element(photoacclimating, state) for state in states(photoacclimating)) ==
        (:carbon, :nitrogen, nothing)

    abundance = Plankton(; states=(abundance=nothing,), reference_state=:abundance)
    @test isnothing(state_element(abundance, reference_state(abundance)))
end

@testset "NiPiZD component declaration" begin
    family = Agate.Models.NiPiZD.NiPiZDFamily()
    components = default_components(family)
    @test keys(components) == (:P, :Z, :N, :D)
    @test components.P isa Plankton
    @test components.Z isa Plankton
    @test all(
        state_element(getproperty(components, name), reference_state(getproperty(components, name))) === :nitrogen
        for name in (:P, :Z)
    )
    @test all(element(getproperty(components, name)) === :nitrogen for name in (:N, :D))
end
