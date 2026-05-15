# # [Size structure example] (@id size_structure_example)

# !!! info
#     This example uses [Oceananigans.jl](https://clima.github.io/OceananigansDocumentation/stable/) and [OceanBioME.jl](https://oceanbiome.github.io/OceanBioME.jl/stable/).
#     We recommend familiarizing yourself with their user interface if you intend to make changes to the physical model setup.

# This example changes the number of plankton and the size structure of the [Agate.jl-NiPiZD](@ref NiPiZD) model. 

# ## Loading dependencies
# The example uses Agate.jl, Oceananigans.jl, and OceanBioME.jl for the ocean simulation.
# CairoMakie is used for plotting.

using Agate
using Agate.Introspection: tracer_names
using Agate.Library.Light
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans.Units
using CairoMakie

nothing #hide

# ## Ecosystem model

# First, we define custom size structures for phytoplankton and zooplankton.
# A diameter range defines the number of size classes (`n`), the minimum and maximum
# equivalent spherical diameters, and the spacing used to generate the classes.

phyto_diameters = (n=3, min_esd=1, max_esd=20, splitting=:log_splitting)
zoo_diameters = (n=3, min_esd=10, max_esd=200, splitting=:log_splitting)
nothing #hide

# The model is constructed as in the quick start, with these diameter definitions passed
# to the NiPiZD constructor. This creates three phytoplankton tracers and three zooplankton
# tracers.

bgc = Agate.Models.NiPiZD.construct(; phyto_diameters, zoo_diameters)
nothing #hide

# We can inspect the tracer names required by the model.

println(tracer_names(bgc))
nothing #hide

# ## Physical model

# Next, we define the light model and combine it with the Agate.jl ecosystem model using OceanBioME.jl.
# The `BoxModelGrid` represents a well-mixed zero-dimensional water column.

light_attenuation = FunctionFieldPAR(; grid=BoxModelGrid())
nothing #hide

bgc_model = Biogeochemistry(bgc; light_attenuation)
nothing #hide

full_model = BoxModel(; biogeochemistry=bgc_model)
nothing #hide

# ## Initial conditions

set!(
    full_model;
    N=8.0,
    P1=0.01,
    P2=0.05,
    P3=0.1,
    Z1=0.01,
    Z2=0.01,
    Z3=0.01,
    D=0.01,
)

# ## Simulation

filename = "size_structure.jld2"

simulation = Simulation(full_model; Δt=240minutes, stop_time=1095days)
nothing #hide

simulation.output_writers[:fields] = JLD2Writer(
    full_model,
    full_model.fields;
    filename,
    schedule=TimeInterval(1day),
    overwrite_existing=true,
)

run!(simulation)
nothing #hide

# ## Plotting

# Load the simulated tracer time series from disk.

tracer_syms = tracer_names(bgc)

timeseries = (; (s => FieldTimeSeries(filename, string(s))[1, 1, 1, :] for s in tracer_syms)...)
nothing #hide

# Plot each tracer concentration through time.

times = FieldTimeSeries(filename, string(first(tracer_syms))).times

fig = Figure(; size=(1200, 1000), fontsize=24)

for (idx, sym) in enumerate(tracer_syms)
    ax = Axis(
        fig[floor(Int, (idx - 1) / 2), Int((idx - 1) % 2)];
        ylabel=string(sym),
        xlabel="Days",
        title="$(sym) concentration (mmol N / m³)",
    )
    lines!(ax, times / day, getproperty(timeseries, sym); linewidth=3)
end

fig