using Agate
using Test

using Agate.Utils: parse_community, class, resolve_class, class_count, build_tracer_index, Tracers
using Agate.Interface: default_plankton_dynamics, default_biogeochem_dynamics, default_community, default_roles

@testset "ClassRef + Tracers accessors" begin
    factory = Agate.Models.NiPiZD.NiPiZDFactory()
    community = default_community(factory)
    plankton_dyn = default_plankton_dynamics(factory)
    biogeochem_dyn = default_biogeochem_dynamics(factory)
    roles = default_roles(factory)

    ctx = parse_community(factory, Float64, community;
                          plankton_dynamics=plankton_dyn,
                          biogeochem_dynamics=biogeochem_dyn,
                          roles=roles)

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
