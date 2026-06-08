using Agate
using Test
using Adapt

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

    active_growth = Agate.Runtime.active_parameters(base_bgc; maximum_growth_rate = (:P1,))
    p = [mu0]
    bgc_p = Agate.Runtime.parameterized(
        base_bgc,
        p;
        active_parameters=active_growth,
    )

    args = (0, 0, 0, 0, 7.0, 1.0, 0.05, 0.05, 0.01, 0.01, 100.0)
    base_tendency = base_bgc(Val(:P1), args...)
    active_tendency = bgc_p(Val(:P1), args...)

    @test active_tendency ≈ base_tendency
    @test bgc_p.parameters.maximum_growth_rate[3] == p[1]
    @test bgc_p.parameters.maximum_growth_rate[4] == base_bgc.parameters.maximum_growth_rate[4]

    active_scalar = Agate.Runtime.active_parameters(base_bgc;
        maximum_growth_rate = (:P1,),
        detritus_remineralization = true,
    )
    p_scalar = [mu0, 2base_bgc.parameters.detritus_remineralization]
    bgc_scalar = Agate.Runtime.parameterized(
        base_bgc,
        p_scalar;
        active_parameters=active_scalar,
    )

    @test bgc_scalar.parameters.maximum_growth_rate[3] == p_scalar[1]
    @test bgc_scalar.parameters.detritus_remineralization == p_scalar[2]
    @test bgc_scalar.parameters.maximum_growth_rate[4] == base_bgc.parameters.maximum_growth_rate[4]


    active_matrix = Agate.Runtime.active_parameters(base_bgc;
        maximum_growth_rate = (:P1,),
        palatability_matrix = ((:Z1, :P1),),
        assimilation_matrix = ((:Z1, :P1),),
    )
    p_matrix = [mu0, 3base_bgc.parameters.interactions.palatability[1, 1], 0.8]
    bgc_matrix = Agate.Runtime.parameterized(
        base_bgc,
        p_matrix;
        active_parameters=active_matrix,
    )

    @test bgc_matrix.parameters.maximum_growth_rate[3] == p_matrix[1]
    @test bgc_matrix.parameters.palatability_matrix[1, 1] == p_matrix[2]
    @test bgc_matrix.parameters.assimilation_matrix[1, 1] == p_matrix[3]
    @test bgc_matrix.parameters.interactions.palatability[1, 1] == p_matrix[2]
    @test bgc_matrix.parameters.interactions.assimilation[1, 1] == p_matrix[3]
    @test bgc_matrix.parameters.interactions.palatability[1, 2] == base_bgc.parameters.interactions.palatability[1, 2]

    matrix_tendency = bgc_matrix(Val(:P1), args...)
    @test matrix_tendency != base_tendency

    p_fast = [2mu0]
    fast_bgc = Agate.Runtime.parameterized(
        base_bgc,
        p_fast;
        active_parameters=active_growth,
    )
    fast_tendency = fast_bgc(Val(:P1), args...)

    @test fast_tendency != base_tendency
    @test fast_bgc.bgc === base_bgc
    @test required_biogeochemical_tracers(fast_bgc) == required_biogeochemical_tracers(base_bgc)
    @test required_biogeochemical_auxiliary_fields(fast_bgc) == required_biogeochemical_auxiliary_fields(base_bgc)
    @test biogeochemical_drift_velocity(fast_bgc, Val(:P1)) == biogeochemical_drift_velocity(base_bgc, Val(:P1))

    adapted_bgc = Adapt.adapt(identity, fast_bgc)
    @test adapted_bgc(Val(:P1), args...) ≈ fast_tendency
end


@testset "ODE problem active parameters" begin
    mu0 = 0.7 / day
    base_bgc = ActiveParameterNiPiZD.construct(;
        parameters=(; maximum_growth_rate=(P1=mu0, P2=mu0)),
    )

    u0 = [7.0, 1.0, 0.05, 0.05, 0.01, 0.01]
    active_growth = Agate.Runtime.active_parameters(base_bgc; maximum_growth_rate = (:P1,))
    p = [mu0]

    problem = Agate.Runtime.ode_problem(
        base_bgc,
        u0,
        (0.0, day);
        p,
        active_parameters=active_growth,
        auxiliary=(; PAR=100.0),
    )

    du = similar(u0)
    problem.f(du, u0, p, 0.0)

    expected = Agate.Runtime.parameterized(
        base_bgc,
        p;
        active_parameters=active_growth,
    )(Val(:P1), 0, 0, 0, 0.0, u0..., 100.0)

    @test du[5] ≈ expected

    p_fast = [2mu0]
    du_fast = similar(u0)
    problem.f(du_fast, u0, p_fast, 0.0)

    @test du_fast[5] != du[5]

    p_scalar = [mu0, 2base_bgc.parameters.detritus_remineralization]
    active_scalar = Agate.Runtime.active_parameters(base_bgc;
        maximum_growth_rate = (:P1,),
        detritus_remineralization = true,
    )
    problem_scalar = Agate.Runtime.ode_problem(
        base_bgc,
        u0,
        (0.0, day);
        p=p_scalar,
        active_parameters=active_scalar,
        auxiliary=(; PAR=100.0),
    )

    du_scalar = similar(u0)
    problem_scalar.f(du_scalar, u0, p_scalar, 0.0)

    expected_scalar = Agate.Runtime.parameterized(
        base_bgc,
        p_scalar;
        active_parameters=active_scalar,
    )(Val(:N), 0, 0, 0, 0.0, u0..., 100.0)

    @test du_scalar[1] ≈ expected_scalar

    active_matrix = Agate.Runtime.active_parameters(base_bgc;
        maximum_growth_rate = (:P1,),
        palatability_matrix = ((:Z1, :P1),),
    )
    p_matrix = [mu0, 3base_bgc.parameters.interactions.palatability[1, 1]]
    problem_matrix = Agate.Runtime.ode_problem(
        base_bgc,
        u0,
        (0.0, day);
        p=p_matrix,
        active_parameters=active_matrix,
        auxiliary=(; PAR=100.0),
    )

    du_matrix = similar(u0)
    problem_matrix.f(du_matrix, u0, p_matrix, 0.0)

    expected_matrix = Agate.Runtime.parameterized(
        base_bgc,
        p_matrix;
        active_parameters=active_matrix,
    )(Val(:P1), 0, 0, 0, 0.0, u0..., 100.0)

    @test du_matrix[5] ≈ expected_matrix
    @test du_matrix[5] != du[5]
end

@testset "active parameter selector validation" begin
    base_bgc = ActiveParameterNiPiZD.construct()

    @test_throws ArgumentError Agate.Runtime.active_parameters(base_bgc; specificity = (:Z1,))
    @test_throws ArgumentError Agate.Runtime.active_parameters(base_bgc; optimum_predator_prey_ratio = (:Z1,))
    @test_throws ArgumentError Agate.Runtime.active_parameters(base_bgc; palatability_matrix = true)
    @test_throws ArgumentError Agate.Runtime.active_parameters(base_bgc; maximum_growth_rate = ((:P1, :P2, :extra),))
    @test_throws ArgumentError Agate.Runtime.active_parameters(base_bgc; detritus_remineralization = false)
    @test_throws ArgumentError Agate.Runtime.active_parameters(base_bgc; not_a_parameter = true)
end


@testset "active parameter set" begin
    base_bgc = ActiveParameterNiPiZD.construct()

    active = Agate.Runtime.active_parameters(base_bgc;
        maximum_growth_rate = (:P1, :P2),
        detritus_remineralization = true,
        palatability_matrix = ((:Z1, :P1), (:Z1, :P2), (:Z2, :P1)),
        assimilation_matrix = ((:Z1, :P1),),
    )

    @test active.layout == (;
        maximum_growth_rate = (; P1 = 1, P2 = 2),
        detritus_remineralization = 3,
        palatability_matrix = (; Z1 = (; P1 = 4, P2 = 5), Z2 = (; P1 = 6)),
        assimilation_matrix = (; Z1 = (; P1 = 7)),
    )

    @test active.labels == (
        "maximum_growth_rate.P1",
        "maximum_growth_rate.P2",
        "detritus_remineralization",
        "palatability_matrix[Z1, P1]",
        "palatability_matrix[Z1, P2]",
        "palatability_matrix[Z2, P1]",
        "assimilation_matrix[Z1, P1]",
    )

    @test length(active) == 7
    @test active.values[1] == base_bgc.parameters.maximum_growth_rate[3]
    @test active.values[2] == base_bgc.parameters.maximum_growth_rate[4]
    @test active.values[3] == base_bgc.parameters.detritus_remineralization
    @test active.values[4] == base_bgc.parameters.interactions.palatability[1, 1]
    @test active.values[5] == base_bgc.parameters.interactions.palatability[1, 2]
    @test active.values[6] == base_bgc.parameters.interactions.palatability[2, 1]
    @test active.values[7] == base_bgc.parameters.interactions.assimilation[1, 1]

    p = copy(active.values)
    p[1] *= 2
    bgc_p = Agate.Runtime.parameterized(base_bgc, p; active_parameters = active)
    @test bgc_p.parameters.maximum_growth_rate[3] == p[1]
end
