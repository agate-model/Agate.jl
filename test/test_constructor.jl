using Agate
using NamedArrays
using Agate.Models.NiPiZD.Tracers
using Agate.Models: NiPiZD

@testset "Models.Constructor" begin
    N2P2ZD_constructed = NiPiZD.construct()

    # declare constants reused in all tests
    P1 = 0.01
    P2 = 0.01
    Z1 = 0.05
    Z2 = 0.05
    N = 7.0
    D = 1
    PAR = 100

    """
    The order in which tracers are called when a Biogeochemistry model is called can vary.
    This function ensures we do not need to hardcode this and so avoid test errors.
    """
    function get_var_order(model)
        # map symbols to variables
        var_map = Dict(:P1 => P1, :P2 => P2, :Z1 => Z1, :Z2 => Z2, :D => D, :N => N)
        # required_biogeochemical_tracers() returns a tuple of symbols
        return [var_map[sym] for sym in required_biogeochemical_tracers(model)]
    end

    @testset "N2P2ZD model" begin

        # N2P2ZD model defined using low level syntax
        include(joinpath("N2P2ZD", "tracers.jl"))
        model = N2P2ZD()
        model_var_order = get_var_order(model)

        # N2P2ZD model constructed from emergent parameters - just using default vals here
        model_constructed = N2P2ZD_constructed()
        model_constructed_var_order = get_var_order(model_constructed)

        @test !iszero(
            model_constructed(Val(:N), 0, 0, 0, 0, model_constructed_var_order..., PAR)
        )
        @test !iszero(
            model_constructed(Val(:D), 0, 0, 0, 0, model_constructed_var_order..., PAR)
        )
        @test !iszero(
            model_constructed(Val(:P1), 0, 0, 0, 0, model_constructed_var_order..., PAR)
        )
        @test !iszero(
            model_constructed(Val(:P2), 0, 0, 0, 0, model_constructed_var_order..., PAR)
        )
        @test !iszero(
            model_constructed(Val(:Z1), 0, 0, 0, 0, model_constructed_var_order..., PAR)
        )
        @test !iszero(
            model_constructed(Val(:Z2), 0, 0, 0, 0, model_constructed_var_order..., PAR)
        )

        @test isapprox(
            model_constructed(Val(:N), 0, 0, 0, 0, model_constructed_var_order..., PAR),
            model(Val(:N), 0, 0, 0, 0, model_var_order..., PAR);
            rtol=0.01,
        )
        @test isapprox(
            model_constructed(Val(:D), 0, 0, 0, 0, model_constructed_var_order..., PAR),
            model(Val(:D), 0, 0, 0, 0, model_var_order..., PAR);
            rtol=0.01,
        )
        @test isapprox(
            model_constructed(Val(:P1), 0, 0, 0, 0, model_constructed_var_order..., PAR),
            model(Val(:P1), 0, 0, 0, 0, model_var_order..., PAR);
            rtol=0.01,
        )
        @test isapprox(
            model_constructed(Val(:P2), 0, 0, 0, 0, model_constructed_var_order..., PAR),
            model(Val(:P2), 0, 0, 0, 0, model_var_order..., PAR);
            rtol=0.01,
        )
        @test isapprox(
            model_constructed(Val(:Z1), 0, 0, 0, 0, model_constructed_var_order..., PAR),
            model(Val(:Z1), 0, 0, 0, 0, model_var_order..., PAR);
            rtol=0.01,
        )
        @test isapprox(
            model_constructed(Val(:Z2), 0, 0, 0, 0, model_constructed_var_order..., PAR),
            model(Val(:Z2), 0, 0, 0, 0, model_var_order..., PAR);
            rtol=0.01,
        )
    end

    @testset "User defined matrices" begin
        names = ["P", "Z"]
        wrong_size_matrix = NamedArray(zeros(Float64, 2, 2), (predator=names, prey=names))
        @test_throws ArgumentError NiPiZD.instantiate(
            N2P2ZD_constructed; palatability_matrix=wrong_size_matrix
        )
        @test_throws ArgumentError NiPiZD.instantiate(
            N2P2ZD_constructed; assimilation_efficiency_matrix=wrong_size_matrix
        )

        # doesn't throw error if dimensions are correct
        names = ["P1", "P2", "Z1", "Z2"]
        correct_size_matrix = NamedArray(zeros(Float64, 4, 4), (predator=names, prey=names))
        new_model = NiPiZD.instantiate(
            N2P2ZD_constructed;
            palatability_matrix=correct_size_matrix,
            assimilation_efficiency_matrix=correct_size_matrix,
        )
    end

    @testset "Diameters passed as an array" begin

        # the default type defined at the top has 2 phyto so expect 2 diameters
        @test_throws ArgumentError NiPiZD.instantiate(
            N2P2ZD_constructed; phyto_diameters=[1]
        )

        # the default type defined at the top has 2 phyto so expect 2 diameters
        @test_throws ArgumentError NiPiZD.instantiate(
            N2P2ZD_constructed; phyto_diameters=[1, 2, 3]
        )

        # diameters can be passed as an array of values rather than a dictionary
        # this is useful in the case where we want 1 phytoplankton with a given diameter
        # it could also be used to fix the diameters of multiple phytoplankton in the model
        NP2ZD = NiPiZD.construct(; n_phyto=1)
        model = NiPiZD.instantiate(NP2ZD; phyto_diameters=[2])
        model_1p_var_order = get_var_order(model)

        # this model only has 1 phyto, 2 zoo tracers (unlike other tests here)
        @test !iszero(model(Val(:N), 0, 0, 0, 0, model_1p_var_order..., PAR))
        @test !iszero(model(Val(:D), 0, 0, 0, 0, model_1p_var_order..., PAR))
        @test !iszero(model(Val(:P1), 0, 0, 0, 0, model_1p_var_order..., PAR))
        @test !iszero(model(Val(:Z1), 0, 0, 0, 0, model_1p_var_order..., PAR))
        @test !iszero(model(Val(:Z2), 0, 0, 0, 0, model_1p_var_order..., PAR))
    end

    @testset "Alternative instantiation" begin

        # N2P2ZD model constructed with user-defined functions (geider growth)
        N2P2ZD_geider = NiPiZD.construct(;
            phyto_args=NiPiZD.DEFAULT_PHYTO_GEIDER_ARGS,
            nutrient_dynamics=nutrients_geider_light,
            phyto_dynamics=phytoplankton_geider_light,
        )

        model_geider = NiPiZD.instantiate(
            N2P2ZD_geider; phyto_args=NiPiZD.DEFAULT_PHYTO_GEIDER_ARGS
        )
        geider_var_order = get_var_order(model_geider)

        @test !iszero(model_geider(Val(:N), 0, 0, 0, 0, geider_var_order..., PAR))
        @test !iszero(model_geider(Val(:D), 0, 0, 0, 0, geider_var_order..., PAR))
        @test !iszero(model_geider(Val(:P1), 0, 0, 0, 0, geider_var_order..., PAR))
        @test !iszero(model_geider(Val(:P2), 0, 0, 0, 0, geider_var_order..., PAR))
        @test !iszero(model_geider(Val(:Z1), 0, 0, 0, 0, geider_var_order..., PAR))
        @test !iszero(model_geider(Val(:Z2), 0, 0, 0, 0, geider_var_order..., PAR))
    end

    @testset "Create objects inside for loop" begin
        prev_p1 = 0
        prev_p2 = 0
        for i in 1:5
            phyto_args = NiPiZD.DEFAULT_PHYTO_ARGS
            phyto_args["allometry"]["maximum_growth_rate"]["a"] = i
            model = NiPiZD.instantiate(N2P2ZD_constructed; phyto_args=phyto_args)
            model_var_order = get_var_order(model)
            p1 = model(Val(:P1), 0, 0, 0, 0, model_var_order..., PAR)
            p2 = model(Val(:P2), 0, 0, 0, 0, model_var_order..., PAR)
            @test !iszero(p1)
            @test !iszero(p2)
            @test model(Val(:P1), 0, 0, 0, 0, model_var_order..., PAR) != prev_p1
            @test model(Val(:P1), 0, 0, 0, 0, model_var_order..., PAR) != prev_p2
            # check the values change as change parameter
            prev_p1 = p1
            prev_p2 = p2
        end
    end

    @testset "Create objects inside function" begin
        function some_wrapper_function(max_growth_rate_a)
            phyto_args = NiPiZD.DEFAULT_PHYTO_ARGS
            phyto_args["allometry"]["maximum_growth_rate"]["a"] = max_growth_rate_a
            model = NiPiZD.instantiate(N2P2ZD_constructed; phyto_args=phyto_args)
            return model
        end

        model1 = some_wrapper_function(5)
        model2 = some_wrapper_function(10)
        model_var_order = get_var_order(model1)

        @test model1(Val(:P1), 0, 0, 0, 0, model_var_order..., PAR) !=
            model2(Val(:P1), 0, 0, 0, 0, model_var_order..., PAR)
        @test model1(Val(:P2), 0, 0, 0, 0, model_var_order..., PAR) !=
            model2(Val(:P2), 0, 0, 0, 0, model_var_order..., PAR)
    end
end
