using Agate
using Test

using Agate.Models: NiPiZDFactory
using Agate.Constructor: update_community, extend_community, update_dynamics, extend_dynamics

@testset "Constructor patch helpers" begin
    factory = NiPiZDFactory()
    community = Agate.FactoryInterface.default_community(factory)
    dynamics = Agate.FactoryInterface.default_plankton_dynamics(factory)

    @testset "update_community (top-level)" begin
        @test_throws ArgumentError update_community(community; X=community.P)

        c2 = update_community(community; P=(n=community.P.n,))
        @test c2.P.n == community.P.n
        @test c2.Z == community.Z

        @test_throws ArgumentError update_community(community; P=(badkey=1,))
        @test_throws ArgumentError update_community(community; P=(pft=community.P.pft,))
    end

    @testset "update_community(group)" begin
        c3 = update_community(community, :P; n=community.P.n)
        @test c3.P.n == community.P.n

        @test_throws ArgumentError update_community(community, :P; badkey=1)
        @test_throws ArgumentError update_community(community, :P; pft=community.P.pft)
    end

    @testset "extend_community" begin
        c4 = extend_community(community; X=community.P)
        @test hasproperty(c4, :X)
        @test c4.X == community.P

        @test_throws ArgumentError extend_community(community; P=community.P)
    end

    @testset "update_dynamics" begin
        d2 = update_dynamics(dynamics; Z=getproperty(dynamics, :Z))
        @test d2 == dynamics

        @test_throws ArgumentError update_dynamics(dynamics; X=1)
    end

    @testset "extend_dynamics" begin
        d3 = extend_dynamics(dynamics; X=getproperty(dynamics, :Z))
        @test hasproperty(d3, :X)

        @test_throws ArgumentError extend_dynamics(dynamics; Z=getproperty(dynamics, :Z))
    end
end
