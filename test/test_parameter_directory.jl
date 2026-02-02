using Agate
using Test

import Agate.Interface: parameter_directory

@testset "Parameter directory" begin
    @testset "NiPiZD" begin
        factory = Agate.Models.NiPiZD.NiPiZDFactory()
        dir = parameter_directory(factory)
        @test !isempty(dir)

        bgc = Agate.NiPiZD.construct(; grid=dummy_grid(Float32))
        dir_names = Set(spec.name for spec in dir)

        # All constructed parameters should be declared in the directory.
        for k in parameter_names(bgc)
            @test k in dir_names
        end

        specmap = Dict(spec.name => spec for spec in dir)
        @test specmap[:detritus_remineralization].shape == :scalar
        @test specmap[:maximum_growth_rate].shape == :vector
        @test specmap[:palatability_matrix].shape == :matrix
        @test specmap[:palatability_matrix].axes == (:consumer, :prey)
        @test specmap[:assimilation_matrix].axes == (:consumer, :prey)
    end

    @testset "DARWIN" begin
        factory = Agate.Models.DARWIN.DarwinFactory()
        dir = parameter_directory(factory)
        @test !isempty(dir)

        bgc = Agate.DARWIN.construct(; grid=dummy_grid(Float32))
        dir_names = Set(spec.name for spec in dir)

        for k in parameter_names(bgc)
            @test k in dir_names
        end

        specmap = Dict(spec.name => spec for spec in dir)
        @test specmap[:DOC_remineralization].shape == :scalar
        @test specmap[:linear_mortality].shape == :vector
        @test specmap[:assimilation_matrix].shape == :matrix
        @test specmap[:palatability_matrix].axes == (:consumer, :prey)
        @test specmap[:assimilation_matrix].axes == (:consumer, :prey)
    end
end
