using Agate
using Agate.Introspection:
    tracer_names,
    parameter_names,
    plankton_groups,
    plankton_tracers,
    nonplankton_tracers,
    tracer_groups
using Test

@testset "Public introspection helpers" begin
    @testset "Model-constructed instance" begin
        bgc = Agate.Models.NiPiZD.construct(; grid=dummy_grid(Float32))

        @test tracer_names(bgc) == [:N, :D, :Z1, :Z2, :P1, :P2]

        groups = tracer_groups(bgc)
        @test groups.all == tracer_names(bgc)
        @test groups.by_group.Z == [:Z1, :Z2]
        @test groups.by_group.P == [:P1, :P2]
        @test groups.plankton == [:Z1, :Z2, :P1, :P2]
        @test groups.nonplankton == [:N, :D]
        @test plankton_groups(bgc) == groups.by_group
        @test plankton_tracers(bgc) == groups.plankton
        @test nonplankton_tracers(bgc) == groups.nonplankton

        pars = parameter_names(bgc)
        @test !isempty(pars)
        @test :data ∉ pars
    end

    @testset "Generated model (define_tracer_functions)" begin
        include(joinpath(@__DIR__, "NPZD", "tracers.jl"))

        model = AgateNPZD(parameters)
        @test tracer_names(model) == [:N, :D, :P, :Z]
        @test plankton_groups(model) == NamedTuple()
        @test plankton_tracers(model) == Symbol[]
        @test nonplankton_tracers(model) == tracer_names(model)

        groups = tracer_groups(model)
        @test groups.all == tracer_names(model)
        @test groups.by_group == NamedTuple()
        @test groups.plankton == Symbol[]
        @test groups.nonplankton == tracer_names(model)
        @test !isempty(parameter_names(model))
    end
end
