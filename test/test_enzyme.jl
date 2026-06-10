using Agate
using Enzyme
using Test

using Oceananigans.Units: day

const EnzymeNiPiZD = Agate.Models.NiPiZD

@testset "Enzyme parameterized tendency gradients" begin
    mu0 = 0.7 / day
    base_bgc = EnzymeNiPiZD.construct(;
        parameters=(; maximum_growth_rate=(P1=mu0, P2=mu0)),
    )

    active = Agate.Runtime.active_parameters(base_bgc;
        maximum_growth_rate = (:P1,),
        detritus_remineralization = true,
        palatability_matrix = ((:Z1, :P1),),
        assimilation_matrix = ((:Z1, :P1),),
    )

    args = (0, 0, 0, 0, 7.0, 1.0, 0.05, 0.05, 0.01, 0.01, 100.0)

    function diagnostic(p)
        parameters = Agate.Runtime.ActiveParameters(base_bgc.parameters, p, active.map)
        p1 = Agate.Runtime.evaluate_tendency(base_bgc, parameters, Val(:P1), args...)
        z1 = Agate.Runtime.evaluate_tendency(base_bgc, parameters, Val(:Z1), args...)
        d = Agate.Runtime.evaluate_tendency(base_bgc, parameters, Val(:D), args...)
        return p1 + 0.5z1 + 0.25d
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
