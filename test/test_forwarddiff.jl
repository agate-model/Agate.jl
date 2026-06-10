using Agate
using Test
using ForwardDiff

const ForwardDiffNiPiZD = Agate.Models.NiPiZD
const ForwardDiffDARWIN = Agate.Models.DARWIN

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

@testset "ForwardDiff NiPiZD ode_problem active parameter smoke test" begin
    mu0 = 0.7 / day
    base_bgc = ForwardDiffNiPiZD.construct(;
        parameters=(; maximum_growth_rate=(P1=mu0, P2=mu0)),
    )

    active_growth = Agate.Runtime.active_parameters(base_bgc; maximum_growth_rate=(:P1,))
    u0 = [7.0, 1.0, 0.05, 0.05, 0.01, 0.01]
    problem = Agate.Runtime.ode_problem(
        base_bgc,
        u0,
        (0.0, day);
        p=copy(active_growth.values),
        active_parameters=active_growth,
        auxiliary=(; PAR=100.0),
    )

    function p1_tendency_with_active_growth_rate(mu)
        T = typeof(mu)
        u = T.(u0)
        du = similar(u)
        problem.f(du, u, [mu], zero(T))
        return du[5]
    end

    dP1_dmu = ForwardDiff.derivative(p1_tendency_with_active_growth_rate, mu0)
    @test isfinite(dP1_dmu)

    δ = mu0 * 1e-6
    fd = (
        p1_tendency_with_active_growth_rate(mu0 + δ) -
        p1_tendency_with_active_growth_rate(mu0 - δ)
    ) / (2δ)
    @test isapprox(dP1_dmu, fd; rtol=1e-4, atol=1e-10)
end

@testset "ForwardDiff DARWIN tendency smoke tests" begin
    function p1_tendency_with_growth_rate(mu)
        T = typeof(mu)
        bgc = ForwardDiffDARWIN.construct(;
            scalar_type=T,
            parameters=(; maximum_growth_rate=T[zero(T), zero(T), mu, T(2.0 / day)]),
        )

        DIC = T(10.0)
        DIN = T(3.0)
        PO4 = T(0.2)
        DOC = T(0.1)
        POC = T(0.1)
        DON = T(0.01)
        PON = T(0.01)
        DOP = T(0.001)
        POP = T(0.001)
        Z1 = T(0.02)
        Z2 = T(0.02)
        P1 = T(0.01)
        P2 = T(0.01)
        PAR = T(100.0)

        return bgc(
            Val(:P1),
            0,
            0,
            0,
            0,
            DIC,
            DIN,
            PO4,
            DOC,
            POC,
            DON,
            PON,
            DOP,
            POP,
            Z1,
            Z2,
            P1,
            P2,
            PAR,
        )
    end

    mu0 = 2.0 / day
    T0 = typeof(mu0)
    bgc0 = ForwardDiffDARWIN.construct(;
        scalar_type=T0,
        parameters=(; maximum_growth_rate=T0[zero(T0), zero(T0), mu0, T0(2.0 / day)]),
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
