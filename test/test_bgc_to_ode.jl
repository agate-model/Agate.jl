using Agate
using Agate.Library.Light

using OrdinaryDiffEq

using Oceananigans.Units

@testset "bgc_to_ode" begin
    include(joinpath("..", "examples", "NPZD", "tracers.jl"))

    init_conditions = (Z=0.05, P=0.01, N=7.0, D=0.0)
    tspan = (0.0, 365day)
    Δt = dt = 7days

    # uses all BGC params by default
    prob = bgc_to_ode(NPZD, cyclical_PAR(; z=-10), init_conditions, tspan)
    @test length(prob.p) == 11

    # can also specify only some BGC parameters - this is useful for inference
    params = (:μ₀, :kₙ, :lᵖᵈ, :α)
    prob2 = bgc_to_ode(NPZD, cyclical_PAR(; z=-10), init_conditions, tspan, params)
    @test length(prob2.p) == 4

    # sanity check solution is correct if only passing some of the BGC parameters
    sol_default = solve(prob, Tsit5(); saveat=Δt)
    # NOTE: using default values in this example as well
    p = [0.6989 / day, 2.3868, 0.0101 / day, 0.1953 / day]
    sol_params = solve(prob2, Tsit5(); p=p, saveat=Δt)

    @test sol_default == sol_params
end
