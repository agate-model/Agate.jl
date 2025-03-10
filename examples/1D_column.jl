# # [1D water column example] (@id 1D_column_example)

# In this example we run a default `Agate.jl-NiPiZD` model inside a simple 1D water column model.
# The physical model setup is based on an example provided in the OceanBioME.jl documentation and represents an idealized 200m deep North Atlantic time series.

# ## Loading dependencies
# The example uses Agate.jl, Oceananigans.jl, and OceanBioME.jl for the ocean simulations.
# CairoMakie is used for plotting.

using Agate
using Agate.Constructors: NiPiZD
using Agate.Library.Light
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans.Units
using CairoMakie

const year = years = 365day
nothing #hide

# ## Ecosystem model

# First, we construct our ecosystem model.
# Here, we use a default 2 phytoplankton, 2 zooplankton `Agate.jl-NiPiZD` ecosystem model.

N2P2ZD = NiPiZD.construct()
nothing #hide

# ## Forcings

# Second, we define the model physical forcings. In this example, mixed layer depth (MLD) forces the physical mixing (diffusivity), while PAR influences plankton photosynthesis.
#diffusivity
@inline function diffusivity(x, y, z, t)
    H(t, t₀, t₁) = ifelse(t₀ < t < t₁, 1.0, 0.0)
    function fmld1(t)
        return H(t, 50days, year) *
               (1 / (1 + exp(-(t - 100days) / 5days))) *
               (1 / (1 + exp((t - 330days) / 25days)))
    end
    function MLD(t)
        return -(
            10 +
            340 * (
                1 - fmld1(year - eps(year)) * exp(-mod(t, year) / 25days) -
                fmld1(mod(t, year))
            )
        )
    end
    return 1e-2 * (1 + tanh((z - MLD(t)) / 10)) / 2 + 1e-4
end

#irradiance
@inline function cyclical_PAR(x, y, z, t)
    PAR⁰ =
        60 *
        (1 - cos((t + 15days) * 2π / year)) *
        (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2
    return PAR⁰ * exp(0.2 * z)
end

#plots
t_range = 0.0:days:(365.0 * days)  # Time range from 0 to 365 days 
z_range = -200.0:10.0:0.0  # Depth range from -200m to 0m 
x, y, z = 0.0, 0.0, 0.0
κₜ_values = [log(diffusivity(x, y, z, t)) for t in t_range, z in z_range]
PAR_values = [log(cyclical_PAR(x, y, z, t)) for t in t_range, z in z_range]

fig_forcing = Figure(; resolution=(1000, 800))
ax1 = Axis(fig_forcing[1, 1]; xlabel="Time (days)", ylabel="PAR0", title="log(irradiance)")
CairoMakie.heatmap!(ax1, t_range ./ days, z_range, PAR_values; colormap=:viridis)

ax2 = Axis(fig_forcing[2, 1]; xlabel="Time (days)", ylabel="Depth (m)", title="log(diffusivity)")
CairoMakie.heatmap!(ax2, t_range ./ days, z_range, κₜ_values; colormap=:viridis)

fig_forcing

# ## Physical model

grid = RectilinearGrid(; size=(1, 1, 25), extent=(20meters, 20meters, 200meters))
nothing #hide

bgc_model = Biogeochemistry(
    N2P2ZD(); light_attenuation=FunctionFieldPAR(; grid, PAR_f=cyclical_PAR)
)
nothing #hide

clock = Clock(; time=0.0)

full_model = NonhydrostaticModel(;
    grid,
    clock,
    closure=ScalarDiffusivity(; ν=diffusivity, κ=diffusivity),
    biogeochemistry=bgc_model,
)
nothing #hide

# ## Initial conditions

set!(full_model; N=7.0, P1=0.01, P2=0.01, Z1=0.05, Z2=0.05, D=0.0)

# ## Simulation
filename = "N2P2ZD_column.jld2"

simulation = Simulation(full_model; Δt=20minutes, stop_time=1years)

simulation.output_writers[:profiles] = JLD2OutputWriter(
    full_model,
    full_model.tracers;
    filename=filename,
    schedule=TimeInterval(1day),
    overwrite_existing=true,
)

run!(simulation)

nothing #hide

# ## Plotting
using CairoMakie

#Load time series data
timeseries = NamedTuple{keys(full_model.tracers)}(
    FieldTimeSeries(filename, "$field") for field in keys(full_model.tracers)
)
nothing #hide

timeseries_keys = keys(timeseries)
nothing #hide

P_keys = filter(k -> startswith(string(k), "P"), timeseries_keys)
Z_keys = filter(k -> startswith(string(k), "Z"), timeseries_keys)

fig = Figure(; size=(1400, 400 * (length(P_keys) + length(Z_keys) + 2)), fontsize=20)

axis_kwargs = (xlabel="Time (days)", ylabel="z (m)", limits=((0, nothing), (-200, 0)))

for (i, key) in enumerate(P_keys)
    x, y, z = nodes(timeseries[key])
    z = collect(z)  # Ensure z is a vector
    times = collect(timeseries[key].times / days)  # Convert times to a vector

    axP = Axis(fig[i, 1]; title="$(key) concentration (mmol N / m³)", axis_kwargs...)
    hmP = heatmap!(
        axP, times, z, Float32.(interior(timeseries[key], 1, 1, :, :)'); colormap=:viridis
    )
    Colorbar(fig[i, 2], hmP)
end

#Plot Zooplankton data
for (i, key) in enumerate(Z_keys)
    x, y, z = nodes(timeseries[key])
    z = collect(z)
    times = collect(timeseries[key].times / days)

    axZ = Axis(
        fig[i + length(P_keys), 1];
        title="$(key) concentration (mmol N / m³)",
        axis_kwargs...,
    )
    hmZ = heatmap!(
        axZ, times, z, Float32.(interior(timeseries[key], 1, 1, :, :)'); colormap=:viridis
    )
    Colorbar(fig[i + length(P_keys), 2], hmZ)
end

#Plot Nutrient data
times = collect(timeseries[:N].times / days)
x, y, z = nodes(timeseries[:N])
z = collect(z)

axD = Axis(
    fig[length(P_keys) + length(Z_keys) + 1, 1];
    title="Nutrient concentration (mmol N / m³)",
    axis_kwargs...,
)
dataD = Float32.(interior(timeseries[:N], 1, 1, :, :)')  # Convert data to Float32
hmD = heatmap!(axD, times, z, dataD; colormap=:viridis)
Colorbar(fig[length(P_keys) + length(Z_keys) + 1, 2], hmD)

#Plot Detritus data
times = collect(timeseries[:D].times / days)
x, y, z = nodes(timeseries[:D])
z = collect(z)

axD = Axis(
    fig[length(P_keys) + length(Z_keys) + 2, 1];
    title="Detritus concentration (mmol N / m³)",
    axis_kwargs...,
)
dataD = Float32.(interior(timeseries[:D], 1, 1, :, :)')  # Convert data to Float32
hmD = heatmap!(axD, times, z, dataD; colormap=:viridis)
Colorbar(fig[length(P_keys) + length(Z_keys) + 2, 2], hmD)

#Save figure
save("N2P2ZD_column.png", fig)

fig  # Display the figure
