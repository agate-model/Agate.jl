using Agate
using Test

using Agate.Constructor: construct, update_community
using Agate.Parameters: parameter_registry, update_registry
using Agate.Utils.Specifications: pft_has
using Agate.Models: NiPiZDFactory, DarwinFactory
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers
using Oceananigans.Units

@testset "Parameters and casting" begin
    @testset "NiPiZD parameter shapes and types" begin
        bgc = construct(NiPiZDFactory(); grid=dummy_grid(Float32))
        p = bgc.parameters

        # Runtime bundle should only include parameters actually referenced by equations.
        @test !hasproperty(p, :diameters)

        # Registry-driven allocation must include the canonical rate parameters.
        @test hasproperty(p, :maximum_growth_rate)
        @test hasproperty(p, :nutrient_half_saturation)
        @test hasproperty(p, :maximum_predation_rate)
        @test hasproperty(p, :detritus_remineralization)

        @test eltype(p.maximum_growth_rate) == Float32
        @test eltype(p.maximum_predation_rate) == Float32
        @test p.detritus_remineralization isa Float32

        # Interaction matrices must exist for default predation dynamics.
        @test hasproperty(p, :palatability_matrix)
        @test hasproperty(p, :assimilation_matrix)
        @test eltype(p.palatability_matrix) == Float32
        @test size(p.palatability_matrix) == (4, 4)
        @test size(p.assimilation_matrix) == (4, 4)

        # Tracer ordering deterministic.
        @test required_biogeochemical_tracers(bgc) == (:N, :D, :Z1, :Z2, :P1, :P2)
    end

    @testset "DARWIN parameter shapes and types" begin
        bgc = construct(DarwinFactory(); grid=dummy_grid(Float32))
        p = bgc.parameters

        # Runtime bundle should not include structural community info.
        @test !hasproperty(p, :diameters)

        # Stoichiometry scalars come from the registry.
        @test hasproperty(p, :nitrogen_to_carbon)
        @test hasproperty(p, :phosphorus_to_carbon)
        @test p.nitrogen_to_carbon isa Float32
        @test p.phosphorus_to_carbon isa Float32

        # Interaction matrices exist for predation terms.
        @test hasproperty(p, :palatability_matrix)
        @test size(p.palatability_matrix) == (4, 4)
    end

    @testset "Registry defaults + overrides" begin
        factory = NiPiZDFactory()
        base = Agate.Models.default_community(factory)

        # Community inputs are structural only: PFTs should not store defaults.
        @test !pft_has(base.Z.pft, :maximum_predation_rate)
        @test !pft_has(base.P.pft, :maximum_growth_rate)

        # Defaults come from registry (zoo predation is non-zero, phyto entries are zero).
        p_default = construct(factory).parameters
        @test p_default.maximum_predation_rate[1] > 0
        @test p_default.maximum_predation_rate[2] > 0
        @test p_default.maximum_predation_rate[3] == 0
        @test p_default.maximum_predation_rate[4] == 0

        # Allometric defaults should vary across sizes (Z1 vs Z2).
        @test p_default.maximum_predation_rate[1] != p_default.maximum_predation_rate[2]

        # Override by updating the registry.
        registry = update_registry(parameter_registry(factory); maximum_predation_rate=(Z=0.5 / day,))
        p_over = construct(factory; grid=dummy_grid(Float32), community=base, registry=registry).parameters

        @test p_over.maximum_predation_rate[1] == Float32(0.5 / day)
        @test p_over.maximum_predation_rate[2] == Float32(0.5 / day)
        @test p_over.maximum_predation_rate[3] == 0f0
        @test p_over.maximum_predation_rate[4] == 0f0
    end

    
@testset "ParamRegistry show" begin
    factory = NiPiZDFactory()
    reg = parameter_registry(factory)

    out = sprint(show, MIME"text/plain"(), reg)
    @test occursin("maximum_growth_rate", out)
    @test occursin("detritus_remineralization", out)
    # Docs should be included in the registry display.
    @test occursin("Maximum phytoplankton growth rate", out) ||
          occursin("Maximum grazing rate", out) ||
          occursin("Detritus remineralization rate", out)

    # Updated registries should also display cleanly.
    reg2 = update_registry(reg; maximum_predation_rate=(Z=0.5 / day,))
    out2 = sprint(show, MIME"text/plain"(), reg2)
    @test occursin("maximum_predation_rate", out2)
end
end
