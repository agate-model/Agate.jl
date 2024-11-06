# ==================================================
# Column example based on
# https://oceanbiome.github.io/OceanBioME.jl/stable/generated/column/
# ==================================================

using Agate
using Agate.Library.Light
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans.Units

const year = years = 365days

# ==================================================
# Define BGC model (NPZD with simple cyclical PAR)
# ==================================================

include(joinpath("NPZD", "tracers.jl"))

# column grid 1x1x40 (depth at increments of 5m)
grid = RectilinearGrid(; size=(1, 1, 40), extent=(20meters, 20meters, 200meters))

bgc_model = Biogeochemistry(
    NPZD(); light_attenuation=FunctionPAR(; grid, PAR_f=cyclical_PAR)
)

# ==================================================
# Specify Oceananigans model
# ==================================================
@inline H(t, t₀, t₁) = ifelse(t₀ < t < t₁, 1.0, 0.0)

@inline fmld1(t) =
    H(t, 50days, year) *
    (1 / (1 + exp(-(t - 100days) / 5days))) *
    (1 / (1 + exp((t - 330days) / 25days)))

@inline MLD(t) = -(
    10 +
    340 * (1 - fmld1(year - eps(year)) * exp(-mod(t, year) / 25days) - fmld1(mod(t, year)))
)

@inline κₜ(x, y, z, t) = 1e-2 * (1 + tanh((z - MLD(t)) / 10)) / 2 + 1e-4

clock = Clock(; time=0.0)

full_model = NonhydrostaticModel(;
    grid, clock, closure=ScalarDiffusivity(; ν=κₜ, κ=κₜ), biogeochemistry=bgc_model
)

# ==================================================
# Set initial conditions
# ==================================================

set!(full_model; N=7.0, P=0.01, Z=0.05, D=0.0)

# ==================================================
# Simulate
# ==================================================

filename = "column.jld2"

simulation = Simulation(full_model; Δt=5minutes, stop_time=3years)
simulation.output_writers[:profiles] = JLD2OutputWriter(
    full_model,
    full_model.tracers;
    filename=filename,
    schedule=TimeInterval(1day),
    overwrite_existing=true,
)

run!(simulation)

timeseries = NamedTuple{keys(full_model.tracers)}(
    FieldTimeSeries(filename, "$field") for field in keys(full_model.tracers)
)

# ==================================================
# Plotting
# ==================================================

using CairoMakie

x, y, z = nodes(timeseries.P)
times = timeseries.P.times

fig = Figure(; size=(1000, 600), fontsize=20)

axis_kwargs = (
    xlabel="Time (days)", ylabel="z (m)", limits=((0, times[end] / days), (-200meters, 0))
)

axP = Axis(fig[1, 1]; title="Phytoplankton concentration (mmol N / m³)", axis_kwargs...)
hmP = CairoMakie.heatmap!(
    times / days, z, interior(timeseries.P, 1, 1, :, :)'; colormap=:batlow
)
Colorbar(fig[1, 2], hmP)

axZ = Axis(fig[2, 1]; title="Zooplankton concentration (mmol N / m³)", axis_kwargs...)
hmZ = CairoMakie.heatmap!(
    times / days, z, interior(timeseries.Z, 1, 1, :, :)'; colormap=:batlow
)
Colorbar(fig[2, 2], hmZ)

fig

save("NPZD_column.png", fig)
