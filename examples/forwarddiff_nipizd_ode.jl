using Agate
using ForwardDiff
using OrdinaryDiffEq: ODEProblem, Tsit5, solve
using CairoMakie
using Printf

using Oceananigans.Biogeochemistry: required_biogeochemical_tracers
using Oceananigans.Units: day

const NiPiZD = Agate.Models.NiPiZD
const TRACERS = (:N, :D, :Z1, :Z2, :P1, :P2)

function nipizd_model_with_p1_growth_rate(mu)
    T = typeof(mu)
    return NiPiZD.construct(;
        scalar_type=T,
        parameters=(; maximum_growth_rate=[zero(T), zero(T), mu, T(0.7 / day)]),
    )
end

function initial_conditions(::Type{T}) where {T}
    return T[7.0, 1.0, 0.05, 0.05, 0.01, 0.01]
end

function constant_par(::Type{T}) where {T}
    return T(100.0)
end

function solve_nipizd(mu; saveat=range(0.0, 30day; length=121))
    T = typeof(mu)
    bgc = nipizd_model_with_p1_growth_rate(mu)

    required_biogeochemical_tracers(bgc) == TRACERS ||
        error("Unexpected NiPiZD tracer order: $(required_biogeochemical_tracers(bgc))")

    function rhs!(du, u, _, t)
        PAR = constant_par(T)
        for (i, tracer) in enumerate(TRACERS)
            du[i] = bgc(Val(tracer), 0, 0, 0, t, u..., PAR)
        end
        return nothing
    end

    u0 = initial_conditions(T)
    problem = ODEProblem(rhs!, u0, (first(saveat), last(saveat)))
    return solve(problem, Tsit5(); saveat=saveat, abstol=1e-10, reltol=1e-10)
end

function phytoplankton_trajectory(theta; saveat=range(0.0, 30day; length=121))
    sol = solve_nipizd(theta[1]; saveat=saveat)
    values = reduce(hcat, sol.u)
    return vec(values[5:6, :])
end

function finite_difference_trajectory(mu0, delta; saveat)
    plus = reshape(phytoplankton_trajectory([mu0 + delta]; saveat=saveat), 2, :)
    minus = reshape(phytoplankton_trajectory([mu0 - delta]; saveat=saveat), 2, :)
    return (plus .- minus) ./ (2delta)
end

saveat = collect(range(0.0, 30day; length=121))
mu0 = 0.7 / day
delta = mu0 * 1e-6

baseline = reshape(phytoplankton_trajectory([mu0]; saveat=saveat), 2, :)
J = ForwardDiff.jacobian(theta -> phytoplankton_trajectory(theta; saveat=saveat), [mu0])
forwarddiff_sensitivity = reshape(J[:, 1], 2, :)
finite_difference_sensitivity = finite_difference_trajectory(mu0, delta; saveat=saveat)

final_forwarddiff = forwarddiff_sensitivity[:, end]
final_finite_difference = finite_difference_sensitivity[:, end]

@printf("Final dP1/dmu, ForwardDiff:       %.8e\n", final_forwarddiff[1])
@printf("Final dP1/dmu, finite difference: %.8e\n", final_finite_difference[1])
@printf("Final dP2/dmu, ForwardDiff:       %.8e\n", final_forwarddiff[2])
@printf("Final dP2/dmu, finite difference: %.8e\n", final_finite_difference[2])
@printf(
    "Maximum absolute sensitivity difference: %.8e\n",
    maximum(abs.(forwarddiff_sensitivity .- finite_difference_sensitivity)),
)

fig = Figure(; size=(900, 620))
time_days = saveat ./ day

ax1 = Axis(fig[1, 1]; xlabel="time (days)", ylabel="phytoplankton concentration", title="NiPiZD phytoplankton trajectory")
lines!(ax1, time_days, baseline[1, :]; label="P1")
lines!(ax1, time_days, baseline[2, :]; label="P2")
axislegend(ax1; position=:rb)

ax2 = Axis(fig[2, 1]; xlabel="time (days)", ylabel="sensitivity to P1 growth rate", title="ForwardDiff sensitivity with finite-difference check")
lines!(ax2, time_days, forwarddiff_sensitivity[1, :]; label="dP1/dmu, ForwardDiff")
lines!(ax2, time_days, finite_difference_sensitivity[1, :]; linestyle=:dash, label="dP1/dmu, finite difference")
lines!(ax2, time_days, forwarddiff_sensitivity[2, :]; label="dP2/dmu, ForwardDiff")
lines!(ax2, time_days, finite_difference_sensitivity[2, :]; linestyle=:dash, label="dP2/dmu, finite difference")
axislegend(ax2; position=:rb)

output_path = joinpath(@__DIR__, "forwarddiff_nipizd_ode_sensitivity.png")
save(output_path, fig)
@info "Saved figure" output_path
