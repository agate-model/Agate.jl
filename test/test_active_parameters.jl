using Agate
using Test

using Oceananigans.Biogeochemistry:
    biogeochemical_drift_velocity,
    required_biogeochemical_auxiliary_fields,
    required_biogeochemical_tracers
using Oceananigans.Units: day

const ActiveParameterNiPiZD = Agate.Models.NiPiZD

@testset "parameterized BGC" begin
    mu0 = 0.7 / day
    base_bgc = ActiveParameterNiPiZD.construct(;
        parameters=(; maximum_growth_rate=(P1=mu0, P2=mu0)),
    )

    p = [mu0]
    bgc_p = Agate.Runtime.parameterized(
        base_bgc,
        p;
        active_parameters=(; maximum_growth_rate=(; P1=1)),
    )

    args = (0, 0, 0, 0, 7.0, 1.0, 0.05, 0.05, 0.01, 0.01, 100.0)
    base_tendency = base_bgc(Val(:P1), args...)
    active_tendency = bgc_p(Val(:P1), args...)

    @test active_tendency ≈ base_tendency
    @test bgc_p.parameters.maximum_growth_rate[3] == p[1]
    @test bgc_p.parameters.maximum_growth_rate[4] == base_bgc.parameters.maximum_growth_rate[4]

    p_fast = [2mu0]
    fast_bgc = Agate.Runtime.parameterized(
        base_bgc,
        p_fast;
        active_parameters=(; maximum_growth_rate=(; P1=1)),
    )
    fast_tendency = fast_bgc(Val(:P1), args...)

    @test fast_tendency != base_tendency
    @test fast_bgc.bgc === base_bgc
    @test required_biogeochemical_tracers(fast_bgc) == required_biogeochemical_tracers(base_bgc)
    @test required_biogeochemical_auxiliary_fields(fast_bgc) == required_biogeochemical_auxiliary_fields(base_bgc)
    @test biogeochemical_drift_velocity(fast_bgc, Val(:P1)) == biogeochemical_drift_velocity(base_bgc, Val(:P1))
end
