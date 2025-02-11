using Agate
using NamedArrays
using Agate.Models.Tracers

@testset "Models.Constructor" begin

    # N2P2ZD model constructed from emergent parameters
    N2P2ZD_constructed = Agate.Constructors.NPZD_size_structured.construct()

    # declare constants reused in all tests
    P1 = 0.01
    P2 = 0.01
    Z1 = 0.05
    Z2 = 0.05
    N = 7.0
    D = 1
    PAR = 100

    @testset "N2P2ZD model" begin
        # N2P2ZD model defined using low level syntax
        include(joinpath("..", "examples", "N2P2ZD", "tracers.jl"))
        model = N2P2ZD()

        model_constructed = Agate.Constructors.NPZD_size_structured.instantiate(
            N2P2ZD_constructed
        )

        @test !iszero(model_constructed(Val(:N), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))
        @test !iszero(model_constructed(Val(:D), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))
        @test !iszero(model_constructed(Val(:P1), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))
        @test !iszero(model_constructed(Val(:P2), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))
        @test !iszero(model_constructed(Val(:Z1), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))
        @test !iszero(model_constructed(Val(:Z2), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))

        @test isapprox(
            model_constructed(Val(:N), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR),
            model(Val(:N), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR);
            rtol=0.01,
        )
        @test isapprox(
            model_constructed(Val(:D), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR),
            model(Val(:D), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR);
            rtol=0.01,
        )
        @test isapprox(
            model_constructed(Val(:P1), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR),
            model(Val(:P1), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR);
            rtol=0.01,
        )
        @test isapprox(
            model_constructed(Val(:P2), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR),
            model(Val(:P2), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR);
            rtol=0.01,
        )
        @test isapprox(
            model_constructed(Val(:Z1), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR),
            model(Val(:Z1), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR);
            rtol=0.01,
        )
        @test isapprox(
            model_constructed(Val(:Z2), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR),
            model(Val(:Z2), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR);
            rtol=0.01,
        )
    end

    @testset "User defined matrices" begin

        names = ["P", "Z"]
        wrong_size_matrix = NamedArray(zeros(Float64, 2, 2), (predator=names, prey=names))
        @test_throws ArgumentError Agate.Constructors.NPZD_size_structured.instantiate(
            N2P2ZD_constructed; palatability_matrix=wrong_size_matrix
        )
        @test_throws ArgumentError Agate.Constructors.NPZD_size_structured.instantiate(
            N2P2ZD_constructed; assimilation_efficiency_matrix=wrong_size_matrix
        )

        # doesn't throw error if dimensions are correct
        names = ["P1", "P2", "Z1", "Z2"]
        correct_size_matrix = NamedArray(zeros(Float64, 4, 4), (predator=names, prey=names))
        new_model = Agate.Constructors.NPZD_size_structured.instantiate(
            N2P2ZD_constructed;
            palatability_matrix=correct_size_matrix,
            assimilation_efficiency_matrix=correct_size_matrix,
        )
    end

    # @testset "Diameters passed as an array" begin

    #     # the default type doe not allow for 1 phytoplankton
    #     # @test_throws ArgumentError Agate.Constructors.NPZD_size_structured.instantiate(
    #     #     N2P2ZD_constructed; phyto_diameters=[1]
    #     #     )

    #     # diameters can be passed as an array of values rather than a dictionary
    #     # this is useful in the case where we want 1 phytoplankton with a given diameter
    #     # it could also be used to fix the diameters of multiple phytoplankton in the model
    #     NP2ZD = Agate.Constructors.NPZD_size_structured.construct(; n_phyto=1)
    #     model = Agate.Constructors.NPZD_size_structured.instantiate(
    #         N2P2ZD; phyto_diameters=[2]
    #         )

    #     # this model only has 1 phyto, 2 zoo tracers (unlike other tests here)
    #     @test !iszero(model(Val(:N), 0, 0, 0, 0, P1, Z1, Z2, N, D, PAR))
    #     @test !iszero(model(Val(:D), 0, 0, 0, 0, P1, Z1, Z2, N, D, PAR))
    #     @test !iszero(model(Val(:P1), 0, 0, 0, 0, P1, Z1, Z2, N, D, PAR))
    #     @test !iszero(model(Val(:Z1), 0, 0, 0, 0, P1, Z1, Z2, N, D, PAR))
    #     @test !iszero(model(Val(:Z2), 0, 0, 0, 0, P1, Z1, Z2, N, D, PAR))


    #     # the default type defined at the top has 2 phyto so expect 2 diameters
    #     @test_throws ArgumentError Agate.Constructors.NPZD_size_structured.instantiate(
    #         N2P2ZD_constructed; phyto_diameters=[1, 2, 3]
    #         )

    # end

    # @testset "Alternative instantiation" begin

    #     # N2P2ZD model constructed with user-defined functions (geider growth)
    #     N2P2ZD_geider = Agate.Constructors.NPZD_size_structured.construct(;
    #         phyto_args=Agate.Constructors.NPZD_size_structured.DEFAULT_PHYTO_GEIDER_ARGS,
    #         nutrient_dynamics=nutrients_geider_light,
    #         phyto_dynamics=phytoplankton_growth_single_nutrient_geider_light,
    #     )

    #     model_geider = Agate.Constructors.NPZD_size_structured.instantiate(N2P2ZD_geider)

    #     @test !iszero(model_geider(Val(:N), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))
    #     @test !iszero(model_geider(Val(:D), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))
    #     @test !iszero(model_geider(Val(:P1), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))
    #     @test !iszero(model_geider(Val(:P2), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))
    #     @test !iszero(model_geider(Val(:Z1), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))
    #     @test !iszero(model_geider(Val(:Z2), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))
    # end

end
