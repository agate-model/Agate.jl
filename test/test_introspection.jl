using Agate
using Test

@testset "Public introspection helpers" begin
    @testset "Model-constructed instance" begin
        bgc = NiPiZD.construct(; grid=dummy_grid(Float32))

        @test tracer_names(bgc) == [:N, :D, :Z1, :Z2, :P1, :P2]

        pars = parameter_names(bgc)
        @test !isempty(pars)
        @test :data ∉ pars
        @test required_parameters(bgc) == pars
    end

    @testset "Generated model (define_tracer_functions)" begin
        include(joinpath(@__DIR__, "NPZD", "tracers.jl"))

        model = AgateNPZD(parameters)
        @test tracer_names(model) == [:N, :D, :P, :Z]
        @test required_parameters(model) == parameter_names(model)
    end
end
