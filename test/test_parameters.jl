using Agate

using Oceananigans.Units

using Agate.Models.NiPiZD.Parameters:
    create_nipizd_parameters

using Agate.Utils:
    DiameterListSpecification, DiameterRangeSpecification


@testset "Models.Parameters" begin
    @testset "FT enforcement and shapes" begin
        p = create_nipizd_parameters(
            Float32;
            n_phyto=2,
            n_zoo=2,
            phyto_diameters=DiameterRangeSpecification(2, 10, :log_splitting),
            zoo_diameters=DiameterRangeSpecification(20, 100, :linear_splitting),
        )

        @test p.n_P == 2
        @test p.n_Z == 2

        @test eltype(p.diameters) == Float32
        @test eltype(p.maximum_growth_rate) == Float32
        @test eltype(p.palatability_matrix) == Float32

        @test length(p.diameters) == 4
        @test length(p.maximum_growth_rate) == 4
        @test size(p.palatability_matrix) == (4, 4)
        @test size(p.assimilation_efficiency_matrix) == (4, 4)
    end

    @testset "Explicit diameter lists" begin
        phyto = DiameterListSpecification([2.0, 5.0])
        zoo = DiameterListSpecification([20.0, 100.0])

        p = create_nipizd_parameters(
            Float32; n_phyto=2, n_zoo=2, phyto_diameters=phyto, zoo_diameters=zoo
        )

        @test p.diameters == Float32[20.0, 100.0, 2.0, 5.0]
    end

    @testset "Units propagate through construction" begin
        p = create_nipizd_parameters(
            Float64;
            n_phyto=1,
            n_zoo=1,
            phyto_diameters=DiameterListSpecification([2.0]),
            zoo_diameters=DiameterListSpecification([50.0]),
        )

        @test isfinite(p.maximum_growth_rate[2])
        @test p.maximum_growth_rate[2] > 0
    end
end
