using Agate
using Test
using ForwardDiff

const ForwardDiffNiPiZD = Agate.Models.NiPiZD

using Oceananigans.Units: day

@testset "ForwardDiff NiPiZD tendency smoke tests" begin
    function p_1_tendency_with_growth_rate(mu)
        T = typeof(mu)
        bgc = ForwardDiffNiPiZD.construct(;
            scalar_type=T,
            parameters=(; maximum_growth_rate=[mu, T(0.7 / day)]),
        )

        N = T(7.0)
        D = T(1.0)
        Z_1 = T(0.05)
        Z_2 = T(0.05)
        P_1 = T(0.01)
        P_2 = T(0.01)
        PAR = T(100.0)

        return bgc(Val(:P_1), 0, 0, 0, 0, N, D, P_1, P_2, Z_1, Z_2, PAR)
    end

    mu0 = 0.7 / day
    T0 = typeof(mu0)
    bgc0 = ForwardDiffNiPiZD.construct(;
        scalar_type=T0,
        parameters=(; maximum_growth_rate=[mu0, T0(0.7 / day)]),
    )
    @test eltype(bgc0.parameters.maximum_growth_rate) === T0

    dP_1_dmu = ForwardDiff.derivative(p_1_tendency_with_growth_rate, mu0)
    @test isfinite(dP_1_dmu)

    δ = mu0 * 1e-6
    fd =
        (p_1_tendency_with_growth_rate(mu0 + δ) - p_1_tendency_with_growth_rate(mu0 - δ)) /
        (2δ)
    @test isapprox(dP_1_dmu, fd; rtol=1e-4, atol=1e-10)
end

@testset "ForwardDiff NiPiZD ode_problem active parameter smoke test" begin
    mu0 = 0.7 / day
    base_bgc = nipizd_growth_fixture(; mu=mu0)

    active_growth = Agate.Runtime.active_parameters(base_bgc; maximum_growth_rate=(:P_1,))
    p1_index = findfirst(==(:P_1), Agate.Introspection.tracer_names(base_bgc))
    u0 = nipizd_u0()
    problem = Agate.Runtime.ode_problem(
        base_bgc,
        u0,
        (0.0, day);
        p=copy(active_growth.values),
        active_parameters=active_growth,
        auxiliary=(; PAR=100.0),
    )

    function p_1_tendency_with_active_growth_rate(mu)
        T = typeof(mu)
        u = T.(u0)
        du = similar(u)
        problem.f(du, u, [mu], zero(T))
        return du[p1_index]
    end

    dP_1_dmu = ForwardDiff.derivative(p_1_tendency_with_active_growth_rate, mu0)
    @test isfinite(dP_1_dmu)

    δ = mu0 * 1e-6
    fd = (
        p_1_tendency_with_active_growth_rate(mu0 + δ) -
        p_1_tendency_with_active_growth_rate(mu0 - δ)
    ) / (2δ)
    @test isapprox(dP_1_dmu, fd; rtol=1e-4, atol=1e-10)
end
