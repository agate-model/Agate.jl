using Agate
using Test

using Agate.Configuration:
    build_plankton_community, parse_community, DiameterRangeSpecification, Population, Pool,
    realize_components
using Agate.Runtime: class, resolve_class, class_count, build_tracer_index, Tracers
using Agate.ModelFamilies: default_components

pool_component_names(family) = Tuple(
    name for (name, component) in pairs(default_components(family)) if
    component isa Agate.Configuration.Pool
)

@testset "Independent interaction roles" begin
    pft = Agate.Configuration.PFTSpecification()
    group = (; diameters=[1.0], pft)
    community = (P=group, B=group, M=group, Z=group)
    context = parse_community(
        Float64,
        community;
        interaction_roles=(consumers=(:Z, :M), prey=(:P, :B, :M)),
    )

    @test context.consumer_indices == [3, 4]
    @test context.prey_indices == [1, 2, 3]
end

@testset "ClassRef + Tracers accessors" begin
    family = Agate.Models.NiPiZD.NiPiZDFamily()
    community = default_nipizd_community()
    ctx = parse_community(
        Float64, community; biogeochem_tracers=pool_component_names(family)
    )

    @test class_count(ctx, :Z) == 2
    @test class_count(ctx, :P) == 2

    @test resolve_class(ctx, class(:Z, 1)) == 1
    @test resolve_class(ctx, class(:Z, 2)) == 2
    @test resolve_class(ctx, class(:P, 1)) == 3
    @test resolve_class(ctx, class(:P, 2)) == 4

    tracer_names = (:N, :D, ctx.plankton_symbols...)
    aux = (:PAR,)
    idx = build_tracer_index(ctx, tracer_names, aux; n_biogeochem_tracers=2)
    tracers = Tracers(idx)

    # A representative kernel-like positional argument tuple.
    args = (10.0, 20.0, 1.0, 2.0, 3.0, 4.0, 42.0)

    @test tracers.N(args) == 10.0
    @test tracers.D(args) == 20.0
    @test tracers.plankton(args, 1) == 1.0
    @test tracers.plankton(args, 2) == 2.0
    @test tracers.plankton(args, 3) == 3.0
    @test tracers.plankton(args, 4) == 4.0
    @test tracers.Z(args, 2) == 2.0
    @test tracers.P(args, 1) == 3.0
    @test tracers.PAR(args) == 42.0

    @test @inferred(tracers.P(args, 2)) == 4.0
end


@testset "Generic component ClassRef" begin
    layout = realize_components((
        B=Population(; currency=:carbon),
        POM=Pool(:carbon; size_structure=[0.5, 5.0, 50.0]),
    ))

    @test class_count(layout, :B) == 1
    @test class_count(layout, :POM) == 3
    @test resolve_class(layout, class(:B, 1)) == 1
    @test resolve_class(layout, class(:POM, 1)) == 2
    @test resolve_class(layout, class(:POM, 3)) == 4
    @test_throws ArgumentError resolve_class(layout, class(:POM, 4))
end

@testset "Diameter input normalization" begin
    family = Agate.Models.NiPiZD.NiPiZDFamily()
    base = default_nipizd_community()

    community = build_plankton_community(
        base;
        diameters=(
            Z=DiameterRangeSpecification(2, 20.0, 100.0, :linear_splitting),
            P=(n=3, min_esd=2.0, max_esd=10.0, splitting=:log_splitting),
        ),
    )

    ctx = parse_community(
        Float64, community; biogeochem_tracers=pool_component_names(family)
    )

    @test class_count(ctx, :Z) == 2
    @test class_count(ctx, :P) == 3
    @test ctx.diameters[1:2] == [20.0, 100.0]
    @test ctx.diameters[3:5] ≈ [2.0, sqrt(20.0), 10.0]
end
