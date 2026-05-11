using Agate
using Test
using ForwardDiff

const ForwardDiffNiPiZD = Agate.Models.NiPiZD

using Oceananigans.Units: day

@testset "ForwardDiff NiPiZD tendency smoke tests" begin
    function p1_tendency_with_growth_rate(mu)
        T = typeof(mu)
        bgc = ForwardDiffNiPiZD.construct(;
            scalar_type=T,
            parameters=(; maximum_growth_rate=[zero(T), zero(T), mu, T(0.7 / day)]),
        )

        N = T(7.0)
        D = T(1.0)
        Z1 = T(0.05)
        Z2 = T(0.05)
        P1 = T(0.01)
        P2 = T(0.01)
        PAR = T(100.0)

        return bgc(Val(:P1), 0, 0, 0, 0, N, D, Z1, Z2, P1, P2, PAR)
    end

    mu0 = 0.7 / day
    T0 = typeof(mu0)
    bgc0 = ForwardDiffNiPiZD.construct(;
        scalar_type=T0,
        parameters=(; maximum_growth_rate=[zero(T0), zero(T0), mu0, T0(0.7 / day)]),
    )
    @test eltype(bgc0.parameters.maximum_growth_rate) === T0

    dP1_dmu = ForwardDiff.derivative(p1_tendency_with_growth_rate, mu0)
    @test isfinite(dP1_dmu)

    δ = mu0 * 1e-6
    fd =
        (p1_tendency_with_growth_rate(mu0 + δ) - p1_tendency_with_growth_rate(mu0 - δ)) /
        (2δ)
    @test isapprox(dP1_dmu, fd; rtol=1e-4, atol=1e-10)
end
