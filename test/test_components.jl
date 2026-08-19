using Test
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers

using Agate.Configuration:
    Population, Pool, currency, size_structure, sinking, realize_components
using Agate.Factories: default_components

@testset "Component authoring" begin
    population = Population(;
        currency=:nitrogen,
        size_structure=(n=3, min_esd=1.0, max_esd=100.0, splitting=:log_splitting),
        sinking=:population_sinking,
    )
    pool = Pool(:carbon; sinking=:pool_sinking)

    @test currency(population) === :nitrogen
    @test size_structure(population).n == 3
    @test sinking(population) === :population_sinking
    @test currency(pool) === :carbon
    @test isnothing(size_structure(pool))
    @test sinking(pool) === :pool_sinking

    layout = realize_components((P=population, D=pool); scalar_type=Float32)
    @test layout.tracer_order == (:P_1, :P_2, :P_3, :D)
    @test layout.component_tracers == (P=(:P_1, :P_2, :P_3), D=(:D,))
    @test layout.component_indices == (P=(1, 2, 3), D=(4,))
    @test collect(layout.component_diameters.P) ≈ Float32[1, 10, 100]
    @test isnothing(layout.component_diameters.D)

    scalar_population = realize_components((B=Population(; currency=:carbon),))
    @test scalar_population.tracer_order == (:B,)
    @test isnothing(scalar_population.component_diameters.B)

    @test_throws ArgumentError Population(; currency=nothing)
    @test_throws ArgumentError Pool(nothing)
    @test_throws ArgumentError realize_components((POM=Pool(:carbon; size_structure=[1, 2]),))
    @test_throws ArgumentError realize_components((P=population, P_1=Pool(:nitrogen)))
end

@testset "NiPiZD component declaration" begin
    factory = Agate.Models.NiPiZD.NiPiZDFactory()
    components = default_components(factory)
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
