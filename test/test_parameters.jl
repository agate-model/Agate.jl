using Agate
using Test

using Agate.Models: construct, NiPiZDFactory, DarwinFactory
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers

@testset "Parameters and casting" begin
    @testset "NiPiZD parameter shapes and eltypes" begin
        bgc = construct(NiPiZDFactory(); FT=Float32)()
        p = bgc.parameters.data

        # Community is 2 zoo + 2 phyto by default.
        @test length(p.diameters) == 4

        # Registry-driven allocation must include the canonical rate parameters.
        @test hasproperty(p, :maximum_growth_rate)
        @test hasproperty(p, :nutrient_half_saturation)
        @test hasproperty(p, :maximum_predation_rate)

        @test eltype(p.diameters) == Float32
        @test eltype(p.maximum_growth_rate) == Float32
        @test eltype(p.maximum_predation_rate) == Float32

        # Interaction matrices must exist for default predation dynamics.
        @test hasproperty(p, :palatability_matrix)
        @test hasproperty(p, :assimilation_efficiency_matrix)
        @test eltype(p.palatability_matrix) == Float32
        @test size(p.palatability_matrix) == (4, 4)
        @test size(p.assimilation_efficiency_matrix) == (4, 4)

        # Tracer ordering deterministic.
        @test required_biogeochemical_tracers(bgc) == (:N, :D, :Z1, :Z2, :P1, :P2)
    end

    @testset "DARWIN parameter shapes and eltypes" begin
        bgc = construct(DarwinFactory(); FT=Float32)()
        p = bgc.parameters.data

        # Default community is 2 zoo + 2 phyto.
        @test length(p.diameters) == 4
        @test eltype(p.diameters) == Float32

        # Stoichiometry scalars come from BiogeochemistrySpecification.
        @test hasproperty(p, :nitrogen_to_carbon)
        @test hasproperty(p, :phosphorus_to_carbon)
        @test p.nitrogen_to_carbon isa Float32
        @test p.phosphorus_to_carbon isa Float32

        # Interaction matrices exist for predation terms.
        @test hasproperty(p, :palatability_matrix)
        @test size(p.palatability_matrix) == (4, 4)
    end
end
