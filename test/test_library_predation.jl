using Test
using NamedArrays
using Agate.Library.Predation  # Assuming the module is in the same file or loaded

# Initialize NamedArrays with both prey and predator names
prey_names = ["P1", "P2"]
predator_names = ["Z1", "Z2"]
all_names = vcat(prey_names, predator_names)

# Initialize NamedArrays
P = NamedArray([10.0, 20.0, 5.0, 8.0], (all_names,))
maximum_predation_rate = NamedArray([0, 0, 1.0, 2.0], (all_names,))
holling_half_saturation = NamedArray([0, 0, 5.0, 10.0], (all_names,))
palatability = NamedArray(
    [
        0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0
        0.9 0.8 0.0 0.0
        0.7 0.6 0.0 0.0
    ],
    (all_names, all_names),
)

assimilation_efficiency = NamedArray(
    [
        0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0
        0.8 0.7 0.0 0.0
        0.6 0.5 0.0 0.0
    ],
    (all_names, all_names),
)

# Nested testset for the entire predation library
@testset "predation library" begin
    # Test holling_type_2
    @testset "holling_type_2" begin
        @test holling_type_2(10.0, 5.0) isa Real
    end

    # Test predation_loss_idealized
    @testset "predation_loss_idealized" begin
        @test predation_loss_idealized(10.0, 5.0, 1.0, 5.0) isa Real
    end

    # Test predation_gain_idealized
    @testset "predation_gain_idealized" begin
        @test predation_gain_idealized(10.0, 5.0, 0.8, 1.0, 5.0) isa Real
    end

    # Test predation_assimilation_loss_idealized
    @testset "predation_assimilation_loss_idealized" begin
        @test predation_assimilation_loss_idealized(10.0, 5.0, 0.8, 1.0, 5.0) isa Real
    end

    # Test predation_loss_preferential
    @testset "predation_loss_preferential" begin
        @test predation_loss_preferential(10.0, 5.0, 1.0, 5.0, 0.9) isa Real
    end

    # Test predation_gain_preferential
    @testset "predation_gain_preferential" begin
        @test predation_gain_preferential(10.0, 5.0, 0.8, 1.0, 5.0, 0.9) isa Real
    end

    # Test predation_assimilation_loss_preferential
    @testset "predation_assimilation_loss_preferential" begin
        @test predation_assimilation_loss_preferential(10.0, 5.0, 0.8, 1.0, 5.0, 0.9) isa
            Real
    end

    # Test summed_predation_loss_preferential
    @testset "summed_predation_loss_preferential" begin
        @test summed_predation_loss_preferential(
            "P1", P, maximum_predation_rate, holling_half_saturation, palatability
        ) isa Real
    end

    # Test summed_predation_gain_preferential
    @testset "summed_predation_gain_preferential" begin
        @test summed_predation_gain_preferential(
            "Z1",
            P,
            assimilation_efficiency,
            maximum_predation_rate,
            holling_half_saturation,
            palatability,
        ) isa Real
    end

    # Test summed_predation_assimilation_loss_preferential
    @testset "summed_predation_assimilation_loss_preferential" begin
        @test summed_predation_assimilation_loss_preferential(
            "Z1",
            P,
            assimilation_efficiency,
            maximum_predation_rate,
            holling_half_saturation,
            palatability,
        ) isa Real
    end

    # Test net_predation_assimilation_loss_preferential
    @testset "net_predation_assimilation_loss_preferential" begin
        @test net_predation_assimilation_loss_preferential(
            P,
            holling_half_saturation,
            maximum_predation_rate,
            assimilation_efficiency,
            palatability,
        ) isa Real
    end

    # Test assimilation_efficiency_emergent_binary
    @testset "assimilation_efficiency_emergent_binary" begin
        prey_data = Dict("can_be_eaten" => 1)
        predator_data = Dict("can_eat" => 1, "assimilation_efficiency" => 0.8)
        @test assimilation_efficiency_emergent_binary(prey_data, predator_data) isa Real
    end
end
