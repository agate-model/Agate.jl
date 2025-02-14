using Test
using NamedArrays
using Agate.Library.Mortality

plankton_names = ["P1", "P2", "Z1", "Z2"]

P = NamedArray([10.0, 20.0, 5.0, 8.0], (plankton_names,))
linear_mortality = NamedArray([0.1, 0.2, 0.05, 0.08], (plankton_names,))
quadratic_mortality = NamedArray([0.01, 0.02, 0.005, 0.008], (plankton_names,))
quota = NamedArray([0.9, 0.8, 0.7, 0.6], (plankton_names,))

DOM_POM_fractionation = 0.5

# Nested testset for the entire mortality module
@testset "mortality module" begin
    @testset "linear_loss" begin
        @test linear_loss(10.0, 0.1) isa Real
    end

    @testset "quadratic_loss" begin
        @test quadratic_loss(10.0, 0.01) isa Real
    end

    @testset "net_linear_loss" begin
        @test net_linear_loss(P, linear_mortality, DOM_POM_fractionation) isa Real
    end

    @testset "net_linear_loss_quota" begin
        @test net_linear_loss_quota(P, linear_mortality, DOM_POM_fractionation, quota) isa
            Real
    end

    @testset "net_quadratic_loss" begin
        @test net_quadratic_loss(P, quadratic_mortality, DOM_POM_fractionation) isa Real
    end

    @testset "net_quadratic_loss_quota" begin
        @test net_quadratic_loss_quota(
            P, quadratic_mortality, DOM_POM_fractionation, quota
        ) isa Real
    end
end
