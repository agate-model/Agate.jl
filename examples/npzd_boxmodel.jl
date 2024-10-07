using Agate

using DifferentialEquations
using Oceananigans
using OceanBioME
using Plots

using Oceananigans.Fields: FunctionField
using Oceananigans.Units

const year = years = 365day

parameters = (
    μ₀ = 0.6989/day,
    kₙ = 2.3868,
    lᵖⁿ = 0.066/day,
    lᶻⁿ = 0.0102/day,
    lᵖᵈ = 0.0101/day,
    gₘₐₓ = 2.1522/day,
    kₚ = 0.5573,
    β = 0.9116,
    lᶻᵈ = 0.3395/day,
    rᵈⁿ = 0.1213/day,
    α = 0.1953/day,
)
aux_field_vars = [:PAR,]

tracers = Dict(
    "N" => :(phytoplankton_metabolic_loss(P, lᵖⁿ)
    + zooplankton_metabolic_loss(Z, lᶻⁿ)
    + remineralization(D, rᵈⁿ)
    - phytoplankton_growth(N, P, PAR, μ₀, kₙ, α)),

    "D" => :(phytoplankton_mortality_loss(P, lᵖᵈ)
    + zooplankton_assimilation_loss(P, Z, β, gₘₐₓ, kₚ)
    + zooplankton_mortality_loss(Z, lᶻᵈ)
    - remineralization(D, rᵈⁿ)),

    "P" => :(phytoplankton_growth(N, P, PAR, μ₀, kₙ, α)
    - phytoplankton_grazing_loss(P, Z, gₘₐₓ, kₚ)
    - phytoplankton_metabolic_loss(P, lᵖⁿ)
    - phytoplankton_mortality_loss(P, lᵖᵈ)),

    "Z" => :(zooplankton_grazing_gain(P, Z, β, gₘₐₓ, kₚ)
    - zooplankton_metabolic_loss(Z, lᶻⁿ)
    - zooplankton_mortality_loss(Z, lᶻᵈ))
)

NPZD = create_bgc_struct(:NPZD, parameters)
add_bgc_methods(NPZD, tracers, auxiliary_fields=aux_field_vars, helper_functions="NPZD_functions.jl")
npzd_model = NPZD()

init_conditions = (N = 7.0, P = 0.01, Z = 0.05, D=0.0)
timeseries = run_boxmodel(npzd_model, init_conditions)

p = plot(timeseries.P, label="P")
plot!(p, timeseries.Z, label="Z")
plot!(p, timeseries.D, label="D")

savefig(p, "NPZD.png")
