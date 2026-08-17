using Agate
using Test

using Agate.Configuration:
    build_plankton_community, parse_community, DiameterRangeSpecification
using Agate.Runtime: class, resolve_class, class_count, build_tracer_index, Tracers
using Agate.Factories: default_biogeochem_dynamics, default_community


struct GenericRoleFixtureFactory <: Agate.Factories.AbstractBGCFactory end

Agate.Construction.recipe_family(::GenericRoleFixtureFactory) = :GenericRoleFixture
Agate.Factories.parameter_definitions(::GenericRoleFixtureFactory) = ()
Agate.Configuration.matrix_definitions(::GenericRoleFixtureFactory) = (;)

@testset "Independent community roles" begin
    pft = Agate.Configuration.PFTSpecification()
    group = (; diameters=[1.0], pft)
    community = (P=group, B=group, M=group, Z=group)
    ecological_roles = (
        phytoplankton=(:P,),
        bacterioplankton=(:B,),
        mixotrophs=(:M,),
        zooplankton=(:Z,),
    )
    interaction_roles = (consumers=(:Z, :M), prey=(:P, :B, :M))
    parameter_roles = (
        producers=(:P, :M), consumers=(:Z, :M), bacterioplankton=(:B,)
    )

    factory = GenericRoleFixtureFactory()
    context = parse_community(
        Float64,
        community;
        interaction_roles,
        parameter_roles,
    )

    @test (
        context.parameter_role_indices,
        context.consumer_indices,
        context.prey_indices,
    ) == (
        (producers=[1, 3], consumers=[3, 4], bacterioplankton=[2]),
        [3, 4],
        [1, 2, 3],
    )

    recipe = Agate.Construction.capture_model_recipe(
        factory;
        community,
        ecological_roles,
        interaction_roles,
        parameter_roles,
        auxiliary_fields=(),
        scalar_type=Float64,
    )
    manifest = Agate.Construction.capture_model_manifest(
        factory,
        (;),
        context;
        tracer_order=Tuple(context.plankton_symbols),
        auxiliary_fields=(),
        ecological_roles,
        explicit_override_keys=(),
        sinking_tracers=nothing,
        open_bottom=true,
        scalar_type=Float64,
    )

    @test (recipe.ecological_roles, recipe.interaction_roles, recipe.parameter_roles) ==
          (ecological_roles, interaction_roles, parameter_roles)
    @test manifest.ecological_roles == ecological_roles
    @test manifest.interaction_role_indices == (consumers=(3, 4), prey=(1, 2, 3))
    @test manifest.parameter_role_indices ==
          (producers=(1, 3), consumers=(3, 4), bacterioplankton=(2,))
end

@testset "ClassRef + Tracers accessors" begin
    factory = Agate.Models.NiPiZD.NiPiZDFactory()
    community = default_community(factory)
    biogeochem_dyn = default_biogeochem_dynamics(factory)
    ctx = parse_community(
        Float64, community; biogeochem_tracers=keys(biogeochem_dyn)
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

@testset "Diameter input normalization" begin
    factory = Agate.Models.NiPiZD.NiPiZDFactory()
    base = default_community(factory)
    biogeochem_dyn = default_biogeochem_dynamics(factory)

    community = build_plankton_community(
        base;
        diameters=(
            Z=DiameterRangeSpecification(2, 20.0, 100.0, :linear_splitting),
            P=(n=3, min_esd=2.0, max_esd=10.0, splitting=:log_splitting),
        ),
    )

    ctx = parse_community(
        Float64, community; biogeochem_tracers=keys(biogeochem_dyn)
    )

    @test class_count(ctx, :Z) == 2
    @test class_count(ctx, :P) == 3
    @test ctx.diameters[1:2] == [20.0, 100.0]
    @test ctx.diameters[3:5] ≈ [2.0, sqrt(20.0), 10.0]
end
