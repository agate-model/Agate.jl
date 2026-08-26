using Agate
using Test
using Adapt

using Agate.Library.Light: FunctionFieldPAR
using OceanBioME: Biogeochemistry, BoxModel, BoxModelGrid

using Oceananigans.Biogeochemistry:
    biogeochemical_drift_velocity,
    required_biogeochemical_auxiliary_fields,
    required_biogeochemical_tracers
using Oceananigans.Units: day

const ActiveParameterNiPiZD = Agate.Models.NiPiZD

@testset "parameterized BGC" begin
    mu0 = 0.7 / day
    base_bgc = ActiveParameterNiPiZD.construct(;
        parameters=(; maximum_growth_rate=(P_1=mu0, P_2=mu0)),
    )

    active_growth = Agate.Runtime.active_parameters(base_bgc; maximum_growth_rate = (:P_1,))
    p = [mu0]
    bgc_p = Agate.Runtime.parameterized(
        base_bgc,
        p;
        active_parameters=active_growth,
    )

    args = (0, 0, 0, 0, 7.0, 1.0, 0.05, 0.05, 0.01, 0.01, 100.0)
    base_tendency = base_bgc(Val(:P_1), args...)
    active_tendency = bgc_p(Val(:P_1), args...)

    @test active_tendency ≈ base_tendency
    @test bgc_p.parameters.maximum_growth_rate[3] == p[1]
    @test bgc_p.parameters.maximum_growth_rate[4] == base_bgc.parameters.maximum_growth_rate[4]

    p_fast = [2mu0]
    fast_bgc = Agate.Runtime.parameterized(
        base_bgc,
        p_fast;
        active_parameters=active_growth,
    )
    fast_tendency = fast_bgc(Val(:P_1), args...)

    @test fast_tendency != base_tendency
    @test fast_bgc.bgc === base_bgc
    @test required_biogeochemical_tracers(fast_bgc) == required_biogeochemical_tracers(base_bgc)
    @test required_biogeochemical_auxiliary_fields(fast_bgc) == required_biogeochemical_auxiliary_fields(base_bgc)
    @test biogeochemical_drift_velocity(fast_bgc, Val(:P_1)) == biogeochemical_drift_velocity(base_bgc, Val(:P_1))

    adapted_bgc = Adapt.adapt(identity, fast_bgc)
    @test adapted_bgc(Val(:P_1), args...) ≈ fast_tendency
end

@testset "parameterized BGC OceanBioME compatibility" begin
    grid = BoxModelGrid()
    mu0 = 0.7 / day
    base_bgc = ActiveParameterNiPiZD.construct(;
        grid,
        parameters=(; maximum_growth_rate=(P_1=mu0, P_2=mu0)),
    )

    active_growth = Agate.Runtime.active_parameters(base_bgc; maximum_growth_rate = (:P_1,))
    p = copy(active_growth.values)
    bgc_p = Agate.Runtime.parameterized(base_bgc, p; active_parameters=active_growth)

    light_attenuation = FunctionFieldPAR(; grid)
    bgc_model = Biogeochemistry(bgc_p; light_attenuation)

    box_model = BoxModel(; biogeochemistry=bgc_model)

    @test !isnothing(box_model)
end


@testset "ODE problem active parameters" begin
    mu0 = 0.7 / day
    base_bgc = ActiveParameterNiPiZD.construct(;
        parameters=(; maximum_growth_rate=(P_1=mu0, P_2=mu0)),
    )

    u0 = [7.0, 1.0, 0.05, 0.05, 0.01, 0.01]
    active_growth = Agate.Runtime.active_parameters(base_bgc; maximum_growth_rate = (:P_1,))
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
    )(Val(:P_1), 0, 0, 0, 0.0, u0..., 100.0)

    @test du[5] ≈ expected

    p_fast = [2mu0]
    du_fast = similar(u0)
    problem.f(du_fast, u0, p_fast, 0.0)

    @test du_fast[5] != du[5]

    p_scalar = [mu0, 2base_bgc.parameters.detritus_remineralization]
    active_scalar = Agate.Runtime.active_parameters(base_bgc;
        maximum_growth_rate = (:P_1,),
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
        maximum_growth_rate = (:P_1,),
        palatability_matrix = ((:Z_1, :P_1),),
    )
    p_matrix = [mu0, 3base_bgc.parameters.palatability_matrix[1, 1]]
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
    )(Val(:P_1), 0, 0, 0, 0.0, u0..., 100.0)

    @test du_matrix[5] ≈ expected_matrix
    @test du_matrix[5] != du[5]
end

function argument_error_message(f)
    err = try
        f()
        nothing
    catch err
        err
    end

    @test err isa ArgumentError
    return sprint(showerror, err)
end

@testset "active parameter selector validation" begin
    base_bgc = ActiveParameterNiPiZD.construct()
    for name in (:specificity, :protection, :optimum_predator_prey_ratio, :assimilation_efficiency)
        @test !hasproperty(base_bgc.parameters, name)
        @test !getproperty(base_bgc.metadata.parameter_axes, name).runtime_bound
    end

    palatability_selectors = (
        () -> Agate.Runtime.active_parameters(base_bgc; specificity = (:Z_1,)),
        () -> Agate.Runtime.active_parameters(base_bgc; protection = (:P_1,)),
        () -> Agate.Runtime.active_parameters(base_bgc; optimum_predator_prey_ratio = (:Z_1,)),
    )

    for selector in palatability_selectors
        @test occursin(":palatability_matrix", argument_error_message(selector))
    end

    assimilation_message = argument_error_message(() ->
        Agate.Runtime.active_parameters(base_bgc; assimilation_efficiency = (:Z_1,))
    )
    @test occursin(":assimilation_matrix", assimilation_message)

    @test_throws ArgumentError Agate.Runtime.active_parameters(base_bgc; palatability_matrix = true)
    @test_throws ArgumentError Agate.Runtime.active_parameters(base_bgc; maximum_growth_rate = ((:P_1, :P_2, :extra),))
    @test_throws ArgumentError Agate.Runtime.active_parameters(base_bgc; detritus_remineralization = false)
    @test_throws ArgumentError Agate.Runtime.active_parameters(base_bgc; not_a_parameter = true)
end


@testset "active parameter set" begin
    base_bgc = ActiveParameterNiPiZD.construct()

    active = Agate.Runtime.active_parameters(base_bgc;
        maximum_growth_rate = (:P_1, :P_2),
        detritus_remineralization = true,
        palatability_matrix = ((:Z_1, :P_1), (:Z_1, :P_2), (:Z_2, :P_1)),
        assimilation_matrix = ((:Z_1, :P_1),),
    )

    @test active.labels == (
        "maximum_growth_rate.P_1",
        "maximum_growth_rate.P_2",
        "detritus_remineralization",
        "palatability_matrix[Z_1, P_1]",
        "palatability_matrix[Z_1, P_2]",
        "palatability_matrix[Z_2, P_1]",
        "assimilation_matrix[Z_1, P_1]",
    )

    @test length(active) == 7
    @test active.values[1] == base_bgc.parameters.maximum_growth_rate[3]
    @test active.values[2] == base_bgc.parameters.maximum_growth_rate[4]
    @test active.values[3] == base_bgc.parameters.detritus_remineralization
    @test active.values[4] == base_bgc.parameters.palatability_matrix[1, 1]
    @test active.values[5] == base_bgc.parameters.palatability_matrix[1, 2]
    @test active.values[6] == base_bgc.parameters.palatability_matrix[2, 1]
    @test active.values[7] == base_bgc.parameters.assimilation_matrix[1, 1]

    p = copy(active.values)
    p[1] *= 2
    p[3] *= 2
    p[4] *= 3
    p[7] = 0.8

    bgc_p = Agate.Runtime.parameterized(base_bgc, p; active_parameters = active)

    @test bgc_p.parameters.maximum_growth_rate[3] == p[1]
    @test bgc_p.parameters.maximum_growth_rate[4] == p[2]
    @test bgc_p.parameters.detritus_remineralization == p[3]
    @test bgc_p.parameters.palatability_matrix[1, 1] == p[4]
    @test bgc_p.parameters.palatability_matrix[1, 2] == p[5]
    @test bgc_p.parameters.palatability_matrix[2, 1] == p[6]
    @test bgc_p.parameters.assimilation_matrix[1, 1] == p[7]

    @test bgc_p.parameters.palatability_matrix[2, 2] == base_bgc.parameters.palatability_matrix[2, 2]
end

@testset "named-group runtime active parameters" begin
    named_bgc = ActiveParameterNiPiZD.construct(;
        size_structure=(;
            phytoplankton=(diat=[2.0], dino=[10.0]),
            zooplankton=(microzoo=[20.0], mesozoo=[100.0]),
        ),
    )
    flat_bgc = ActiveParameterNiPiZD.construct(;
        size_structure=(;
            phytoplankton=(P=[2.0, 10.0],), zooplankton=(Z=[20.0, 100.0],)
        )
    )

    function runtime_active(bgc, producers, consumers)
        producer_1, producer_2 = producers
        consumer_1, consumer_2 = consumers
        return Agate.Runtime.active_parameters(bgc;
            maximum_growth_rate=(producer_1, producer_2),
            palatability_matrix=((consumer_1, producer_1),),
            assimilation_matrix=((consumer_2, producer_2),),
        )
    end

    named_active = runtime_active(
        named_bgc, (:diat_1, :dino_1), (:microzoo_1, :mesozoo_1)
    )
    flat_active = runtime_active(flat_bgc, (:P_1, :P_2), (:Z_1, :Z_2))

    @test named_active.labels == (
        "maximum_growth_rate.diat_1",
        "maximum_growth_rate.dino_1",
        "palatability_matrix[microzoo_1, diat_1]",
        "assimilation_matrix[mesozoo_1, dino_1]",
    )
    @test named_active.values == flat_active.values

    p = copy(named_active.values)
    p[1] *= 1.2
    p[3] *= 0.8
    p[4] *= 0.9

    args = (0, 0, 0, 0, 7.0, 1.0, 0.05, 0.05, 0.01, 0.01, 100.0)
    named_parameterized = Agate.Runtime.parameterized(
        named_bgc, p; active_parameters=named_active
    )
    flat_parameterized = Agate.Runtime.parameterized(
        flat_bgc, p; active_parameters=flat_active
    )
    @test named_parameterized(Val(:diat_1), args...) ≈
          flat_parameterized(Val(:P_1), args...)

    u0 = [7.0, 1.0, 0.05, 0.05, 0.01, 0.01]
    function rhs(bgc, active)
        problem = Agate.Runtime.ode_problem(
            bgc,
            u0,
            (0.0, day);
            p,
            active_parameters=active,
            auxiliary=(; PAR=100.0),
        )
        du = similar(u0)
        problem.f(du, u0, p, 0.0)
        return du
    end

    @test rhs(named_bgc, named_active) ≈ rhs(flat_bgc, flat_active)
end
