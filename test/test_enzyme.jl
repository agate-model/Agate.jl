using Agate
using Enzyme
using Test

using Oceananigans.Units: day

const EnzymeNiPiZD = Agate.Models.NiPiZD

@testset "Enzyme parameterized tendency smoke test" begin
    mu0 = 0.7 / day
    base_bgc = EnzymeNiPiZD.construct(;
        parameters=(; maximum_growth_rate=(P1=mu0, P2=mu0)),
    )

    args = (0, 0, 0, 0, 7.0, 1.0, 0.05, 0.05, 0.01, 0.01, 100.0)
    active = Agate.Runtime.active_parameters(base_bgc; maximum_growth_rate = (:P1,))

    function p1_tendency(p)
        parameters = Agate.Runtime.ActiveParameters(base_bgc.parameters, p, active.map)
        return Agate.Runtime.evaluate_tendency(base_bgc, parameters, Val(:P1), args...)
    end

    p0 = [mu0]
    grad = zeros(length(p0))
    Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse),
                    Enzyme.Const(p1_tendency),
                    Enzyme.Active,
                    Enzyme.Duplicated(p0, grad))

    @test length(grad) == 1
    @test isfinite(grad[1])

    δ = mu0 * 1e-6
    fd = (p1_tendency([mu0 + δ]) - p1_tendency([mu0 - δ])) / (2δ)
    @test isapprox(grad[1], fd; rtol=1e-4, atol=1e-10)
end
