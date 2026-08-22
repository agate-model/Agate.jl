using Agate
using Enzyme
using Test

using Oceananigans.Units: day

const EnzymeNiPiZD = Agate.Models.NiPiZD

@testset "Enzyme parameterized tendency gradients" begin
    mu0 = 0.7 / day
    base_bgc = EnzymeNiPiZD.construct(;
        parameters=(; maximum_growth_rate=(P_1=mu0, P_2=mu0)),
    )

    active = Agate.Runtime.active_parameters(base_bgc;
        maximum_growth_rate = (:P_1,),
        detritus_remineralization = true,
        palatability_matrix = ((:Z_1, :P_1),),
        assimilation_matrix = ((:Z_1, :P_1),),
    )

    args = (0, 0, 0, 0, 7.0, 1.0, 0.05, 0.05, 0.01, 0.01, 100.0)

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
