using Agate
using Agate.Library.Light
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans.Units
using Agate.Constructor: construct_factory
using Agate.Models: NiPiZDFactory
using Agate.Library.Photosynthesis
using Oceananigans.Fields: FunctionField, ConstantField
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers
using CairoMakie

# Generate a CPU instance (used here to query tracer names).
bgc_cpu = construct_factory(NiPiZDFactory())

# IMPORTANT: get tracer names from the CPU instance (Agate defines the method on this type)
tracer_names = required_biogeochemical_tracers(bgc_cpu)

# Construct a GPU-ready instance for embedding in Oceananigans / OceanBioME models.
bgc_instance = construct_factory(NiPiZDFactory(); arch=GPU())

const year = years = 365days
nothing #hide

@inline PAR⁰(x, y, t) =
    60 *
    (1 - cos((t + 15days) * 2π / year)) *
    (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2

@inline PAR_f(x, y, z, t) = PAR⁰(x, y, t) * exp(0.2 * z)

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

@inline temp(x, y, z, t) = 2.4 * cos(t * 2π / year + 50days) + 10

grid = RectilinearGrid(GPU(); size=(1, 1, 50), extent=(20meters, 20meters, 200meters))

# Build light attenuation on the column grid
light_attenuation = FunctionFieldPAR(; grid, PAR_f=PAR_f)

biogeochemistry = Biogeochemistry(bgc_instance; light_attenuation)

clock = Clock(; time=0.0)
T = FunctionField{Center,Center,Center}(temp, grid; clock)
S = ConstantField(35)

model = NonhydrostaticModel(;
    grid,
    clock,
    tracers=tracer_names,  # IMPORTANT: ensure P1/P2/Z1/Z2/N/D exist
    closure=ScalarDiffusivity(; ν=κₜ, κ=κₜ),
    biogeochemistry,
    auxiliary_fields=(; T, S),
)

set!(model; P1=0.01, P2=0.01, Z1=0.05, Z2=0.05, N=7.0, D=1)

simulation = Simulation(model; Δt=3minutes, stop_time=365days)

filename = "column"
simulation.output_writers[:profiles] = JLD2Writer(
    model,
    model.tracers;
    filename="$filename.jld2",
    schedule=TimeInterval(1day),
    overwrite_existing=true,
)

run!(simulation)

#Load time series data
timeseries = NamedTuple{keys(model.tracers)}(
    FieldTimeSeries(filename, "$field") for field in keys(model.tracers)
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
fig = Figure(; size=(800, 1200), fontsize=16)

#Plot all fields
for (i, key) in enumerate(all_keys)
    x_nodes, y_nodes, z_nodes = nodes(timeseries[key])
    z_vals = collect(z_nodes)
    times = collect(timeseries[key].times / days)

    ax = Axis(
        fig[i, 1];
        title="$(key) concentration (mmol N / m³)",
        xlabel="Time (days)",
        ylabel="z (m)",
        limits=((0, 365), (-200, 0)),
    )
    hm = heatmap!(
        ax,
        times,
        z_vals,
        Float32.(interior(timeseries[key], 1, 1, :, :)');
        colormap=:viridis,
        rasterize=true,
    )  # Rasterize for smaller output
    Colorbar(fig[i, 2], hm)
end

#Save figure
save("N2P2ZD_column_GPU.png", fig)

fig  # Display the figure
