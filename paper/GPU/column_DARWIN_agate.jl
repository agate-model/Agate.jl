using Agate
using Agate.Library.Light
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans.Units
using Agate.Models: construct, DarwinFactory
using Agate.Library.Photosynthesis
using Oceananigans, Printf
using Oceananigans.Fields: FunctionField, ConstantField
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers
using Adapt
using CUDA
using CairoMakie

# Generate the biogeochemical model
bgc_type = construct(DarwinFactory(); FT=Float64)

# Create an instance of the model (CPU)
bgc_instance = bgc_type()

# IMPORTANT: get tracer names from the CPU instance (Agate defines the method on this type)
tracer_names = required_biogeochemical_tracers(bgc_instance)

# Test GPU compatibility
adapted_instance = Adapt.adapt(CuArray, bgc_instance)

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

biogeochemistry = Biogeochemistry(adapted_instance; light_attenuation)

clock = Clock(; time=0.0)
T = FunctionField{Center,Center,Center}(temp, grid; clock)
S = ConstantField(35)

model = NonhydrostaticModel(;
    grid,
    clock,
    tracers=tracer_names, # IMPORTANT: ensure all DARWIN tracers exist
    closure=ScalarDiffusivity(; ν=κₜ, κ=κₜ),
    biogeochemistry,
    auxiliary_fields=(; T, S),
)

# A simple initial condition set: carbon pools in mmol C/m³, nitrogen pools in mmol N/m³,
# phosphate pools in mmol P/m³ (units are treated consistently by the model's stoichiometry).
set!(
    model;
    DIC=2200.0,
    DIN=7.0,
    PO4=0.5,
    DOC=0.0,
    POC=0.0,
    DON=0.0,
    PON=0.0,
    DOP=0.0,
    POP=0.0,
    P1=0.01,
    P2=0.01,
    Z1=0.05,
    Z2=0.05,
)

simulation = Simulation(model; Δt=3minutes, stop_time=365days)

filename = "column_DARWIN"
simulation.output_writers[:profiles] = JLD2Writer(
    model,
    model.tracers;
    filename="$filename.jld2",
    schedule=TimeInterval(1day),
    overwrite_existing=true,
)

run!(simulation)

# Load time series data
timeseries = NamedTuple{keys(model.tracers)}(
    FieldTimeSeries(filename, "$field") for field in keys(model.tracers)
)

timeseries_keys = keys(timeseries)
nothing #hide

# Filter keys for plankton fields
P_keys = filter(k -> startswith(string(k), "P"), timeseries_keys)
Z_keys = filter(k -> startswith(string(k), "Z"), timeseries_keys)

# Order remaining tracers in a consistent way
core_keys = [:DIC, :DIN, :PO4, :DOC, :POC, :DON, :PON, :DOP, :POP]

# Combine all keys into a single list for iteration
all_keys = [P_keys..., Z_keys..., core_keys...]

# Create figure with appropriate size
fig = Figure(; size=(800, 220 * length(all_keys)), fontsize=16)

# Plot all fields
for (i, key) in enumerate(all_keys)
    x_nodes, y_nodes, z_nodes = nodes(timeseries[key])
    z_vals = collect(z_nodes)
    times = collect(timeseries[key].times / days)

    ax = Axis(
        fig[i, 1];
        title="$(key) concentration (mmol / m³)",
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
    )
    Colorbar(fig[i, 2], hm)
end

# Save figure
save("DARWIN_column_GPU.png", fig)

fig # Display the figure
