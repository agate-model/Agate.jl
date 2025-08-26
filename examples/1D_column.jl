# # [1D water column example] (@id 1D_column_example)

# !!! info
#     This example uses [Oceananigans.jl](https://clima.github.io/OceananigansDocumentation/stable/) and [OceanBioME.jl](https://oceanbiome.github.io/OceanBioME.jl/stable/).
#     We recommend familiarizing yourself with their user interface if you intend to make changes to the physical model setup.

# In this example we run a default [Agate.jl-NiPiZD](@ref NiPiZD) model inside a simple 1D water column model.
# The physical model setup is based on an example provided in the OceanBioME.jl documentation and represents an idealized 200m deep North Atlantic time series.

# ## Loading dependencies
# The example uses Agate.jl, Oceananigans.jl, and OceanBioME.jl for the ocean simulations.
# CairoMakie is used for plotting.

using Agate
using Agate.Models: NiPiZD
using Agate.Library.Light
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans.Units
using CairoMakie

import Agate.Library.Light: cyclical_PAR

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
κₜ_values = [diffusivity(x, y, z, t) for t in t_range, z in z_range]
PAR_values = [cyclical_PAR(x, y, z, t) for t in t_range, z in z_range]

fig_forcing = Figure(; resolution=(1000, 800))
ax1 = Axis(fig_forcing[1, 1]; xlabel="Time (days)", ylabel="Depth (m)", title="irradiance")
CairoMakie.heatmap!(ax1, t_range ./ days, z_range, PAR_values; colormap=:viridis)

ax2 = Axis(fig_forcing[2, 1]; xlabel="Time (days)", ylabel="Depth (m)", title="diffusivity")
CairoMakie.heatmap!(ax2, t_range ./ days, z_range, κₜ_values; colormap=:viridis)

fig_forcing

# ## Physical model

grid = RectilinearGrid(; size=(1, 1, 25), extent=(20meters, 20meters, 200meters))
nothing #hide

bgc_model = Biogeochemistry(
    N2P2ZD(); light_attenuation=FunctionFieldPAR(; grid, PAR_f=cyclical_PAR)
)
nothing #hide

full_model = NonhydrostaticModel(;
    grid,
    clock=Clock(; time=0.0),
    closure=ScalarDiffusivity(; ν=diffusivity, κ=diffusivity),
    biogeochemistry=bgc_model,
)
nothing #hide

# ## Initial conditions

set!(full_model; N=7.0, P1=0.01, P2=0.01, Z1=0.05, Z2=0.05, D=0.0) # mmol N / m³

# ## Simulation
filename = "N2P2ZD_column.jld2"

simulation = Simulation(full_model; Δt=20minutes, stop_time=1year)

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

#Load time series data
timeseries = NamedTuple{keys(full_model.tracers)}(
    FieldTimeSeries(filename, "$field") for field in keys(full_model.tracers)
)

timeseries_keys = keys(timeseries)
nothing #hide

#Filter keys for P, Z, N, and D fields
P_keys = filter(k -> startswith(string(k), "P"), timeseries_keys)
Z_keys = filter(k -> startswith(string(k), "Z"), timeseries_keys)
N_key = :N
D_key = :D

#Combine all keys into a single list for iteration
all_keys = [P_keys..., Z_keys..., N_key, D_key]

#Create figure with appropriate size
fig = Figure(; size=(250 * length(all_keys), 400 * length(all_keys)), fontsize=20)

#Axis configuration
axis_kwargs = (xlabel="Time (days)", ylabel="z (m)", limits=((0, nothing), (-200, 0)))

#Plot all fields
for (i, key) in enumerate(all_keys)
    x, y, z = nodes(timeseries[key])
    z = collect(z)  # Ensure z is a vector
    times = collect(timeseries[key].times / days)  # Convert times to a vector

    #Create axis and heatmap
    ax = Axis(fig[i, 1]; title="$(key) concentration (mmol N / m³)", axis_kwargs...)
    hm = heatmap!(
        ax, times, z, Float32.(interior(timeseries[key], 1, 1, :, :)'); colormap=:viridis
    )
    Colorbar(fig[i, 2], hm)
end

#Save figure
save("N2P2ZD_column.png", fig)

fig  # Display the figure
