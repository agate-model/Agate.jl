using Agate
using Agate.Introspection: parameter_names
using Test

import Agate.Factories: FillDefault, parameter_definitions, parameter_directory

@testset "Parameter directory" begin
    @testset "NiPiZD" begin
        factory = Agate.Models.NiPiZD.NiPiZDFactory()
        dir = parameter_directory(factory)
        @test !isempty(dir)

        bgc = Agate.Models.NiPiZD.construct(; grid=dummy_grid(Float32))
        dir_names = Set(spec.name for spec in dir)

        # All constructed parameters should be declared in the directory.
        for k in parameter_names(bgc)
            @test k in dir_names
        end

        specmap = Dict(spec.name => spec for spec in dir)
        definitions = Dict(def.spec.name => def for def in parameter_definitions(factory))
        @test definitions[:linear_mortality].default isa FillDefault
        @test all(!isempty(spec.provides) for spec in values(specmap))
        @test length(specmap[:linear_mortality].provides) == 2
        @test length(specmap[:mortality_export_fraction].provides) == 3
        @test specmap[:detritus_remineralization].shape == :scalar
        @test specmap[:maximum_growth_rate].shape == :vector
        @test specmap[:maximum_growth_rate].axes == :plankton
        @test !isnothing(specmap[:maximum_growth_rate].materialization)
        @test specmap[:maximum_growth_rate].materialization.role === :producers
        @test specmap[:maximum_growth_rate].materialization.fill_value == 0
        @test !isnothing(specmap[:linear_mortality].materialization)
        @test isnothing(specmap[:linear_mortality].materialization.role)
        @test specmap[:linear_mortality].materialization.fill_value == 0
        @test specmap[:palatability_matrix].shape == :matrix
        @test specmap[:palatability_matrix].axes == (:consumer, :prey)
        @test specmap[:assimilation_matrix].axes == (:consumer, :prey)
    end
end
