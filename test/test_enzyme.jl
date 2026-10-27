using Agate
using Enzyme
using ForwardDiff
using Test

using Oceananigans.Units: day
using Oceananigans.Biogeochemistry:
    required_biogeochemical_tracers, required_biogeochemical_auxiliary_fields


@testset "Enzyme parameterized tendency gradients" begin
    mu0 = 0.7 / day
    base_bgc = runtime_npz_growth_fixture(; mu=mu0)

    active = Agate.Runtime.active_parameters(base_bgc;
        maximum_growth_rate = (:P_1,),
        detritus_remineralization = true,
        palatability_matrix = ((:Z_1, :P_1),),
        assimilation_matrix = ((:Z_1, :P_1),),
    )

    args = runtime_npz_args()

    function diagnostic(p)
        bgc_p = Agate.Runtime.parameterized(base_bgc, p; active_parameters=active)
        p_1 = bgc_p(Val(:P_1), args...)
        z_1 = bgc_p(Val(:Z_1), args...)
        d = bgc_p(Val(:D), args...)
        return p_1 + 0.5z_1 + 0.25d
    end

    p0 = copy(active.values)
    grad = zeros(length(p0))
    Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse),
                    Enzyme.Const(diagnostic),
                    Enzyme.Active,
                    Enzyme.Duplicated(p0, grad))

    @test length(grad) == length(active)
    @test all(isfinite, grad)

    for i in eachindex(p0)
        δ = max(abs(p0[i]), 1.0) * 1e-6
        p_plus = copy(p0)
        p_minus = copy(p0)
        p_plus[i] += δ
        p_minus[i] -= δ
        fd = (diagnostic(p_plus) - diagnostic(p_minus)) / (2δ)
        @test isapprox(grad[i], fd; rtol=1e-4, atol=1e-10)
    end
end

@testset "Enzyme quota uptake tendency gradient" begin
    bgc = Agate.Construction.construct(quota_definition())
    active = Agate.Runtime.active_parameters(bgc; maximum_nitrogen_uptake=(:P_1,))
    args = quota_runtime_args()
    diagnostic(p) = Agate.Runtime.parameterized(
        bgc, p; active_parameters=active
    )(Val(:P_1_nitrogen), args...)
    gradient = zeros(length(active.values))
    Enzyme.autodiff(
        Enzyme.set_runtime_activity(Enzyme.Reverse), Enzyme.Const(diagnostic), Enzyme.Active,
        Enzyme.Duplicated(copy(active.values), gradient),
    )
    @test length(gradient) == 1 && isfinite(only(gradient)) && !iszero(only(gradient))
end

@testset "Enzyme FrankenLOBSTER Agate-core gradients" begin
    bgc = Agate.Models.FrankenLOBSTER.construct(; scalar_type=Float64).runtime
    active = Agate.Runtime.active_parameters(
        bgc; maximum_growth_rate=(:P_1,), nitrate_ammonia_inhibition=true
    )
    state = (
        NO₃=2.0, NH₄=0.2, DOM=0.1,
        solid_waste=0.0, inorganic_waste=0.0, dissolved_waste=0.0,
        P_1=0.02, P_2=0.01, Z_1=0.005, Z_2=0.005, H_1=0.01,
    )
    drivers = (; PAR=100.0, T=20.0)
    args = (
        0.0, 0.0, 0.0, 0.0,
        (getproperty(state, name) for name in required_biogeochemical_tracers(bgc))...,
        (getproperty(drivers, name) for name in required_biogeochemical_auxiliary_fields(bgc))...,
    )
    diagnostic(p) = Agate.Runtime.parameterized(
        bgc, p; active_parameters=active
    )(Val(:P_1), args...)

    p0 = copy(active.values)
    gradient = zeros(length(p0))
    Enzyme.autodiff(
        Enzyme.set_runtime_activity(Enzyme.Reverse), Enzyme.Const(diagnostic),
        Enzyme.Active, Enzyme.Duplicated(p0, gradient),
    )
    @test all(isfinite, gradient)
    @test gradient[1] > 0 && gradient[2] < 0
    @test isapprox(gradient, ForwardDiff.gradient(diagnostic, p0); rtol=1e-5, atol=1e-10)
end
