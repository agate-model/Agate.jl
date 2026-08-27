using Agate
using Enzyme
using Test

using Oceananigans.Units: day


@testset "Enzyme parameterized tendency gradients" begin
    mu0 = 0.7 / day
    base_bgc = nipizd_growth_fixture(; mu=mu0)

    active = Agate.Runtime.active_parameters(base_bgc;
        maximum_growth_rate = (:P_1,),
        detritus_remineralization = true,
        palatability_matrix = ((:Z_1, :P_1),),
        assimilation_matrix = ((:Z_1, :P_1),),
    )

    args = nipizd_runtime_args()

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
