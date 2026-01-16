# Quick start

!!! info
    
    Agate.jl is designed to interface with [Oceananigans.jl](https://clima.github.io/OceananigansDocumentation/stable/) and [OceanBioME.jl](https://oceanbiome.github.io/OceanBioME.jl/stable/). We thus recommend familiarizing yourself with their user interface.

### Loading dependencies

```@example quickstart

using Agate
using Agate.Library.Light
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans.Units
using CairoMakie
```

### Model setup

First, we construct our ecosystem model.
Here, we use a default 2 phytoplankton, 2 zooplankton `Agate.jl-NiPiZD` ecosystem model.

```@example quickstart

bgc = construct(NiPiZDFactory())
nothing #hide

# Inspect the required tracer names.
println(tracer_names(bgc))
nothing #hide
```

Next, we define a light function, here we use a default seasonal PAR curve:

```@example quickstart
light_attenuation = FunctionFieldPAR(; grid=BoxModelGrid())
nothing #hide
```

These two models are then combined using OceanBioME.jl

```@example quickstart

bgc_model = Biogeochemistry(bgc; light_attenuation=light_attenuation)
nothing #hide

full_model = BoxModel(; biogeochemistry=bgc_model)
nothing #hide
```

And finally simulated using Oceananigans.jl

```@example quickstart

set!(full_model; N=7.0, P1=0.01, Z1=0.01, P2=0.1, Z2=0.01, D=0.01)

simulation = Simulation(full_model; Δt=240minutes, stop_time=1095days)
nothing #hide

filename = "quick_start.jld2"

simulation.output_writers[:fields] = JLD2Writer(
    full_model,
    full_model.fields;
    filename=filename,
    schedule=TimeInterval(1day),
    overwrite_existing=true,
)

run!(simulation)
nothing #hide
```

### Plotting

```@example quickstart

tracer_syms = tracer_names(bgc)

# Extract data for plotting
timeseries = (; (s => FieldTimeSeries(filename, string(s))[1, 1, 1, :] for s in tracer_syms)...)
nothing #hide

# Create a figure
times = FieldTimeSeries(filename, string(first(tracer_syms))).times

fig = Figure(; size=(1200, 800), fontsize=24)

axs = []
for (idx, sym) in enumerate(tracer_syms)
    ax = Axis(
        fig[floor(Int, (idx - 1) / 2), Int((idx - 1) % 2)];
        ylabel=string(sym),
        xlabel="Days",
        title="$(sym) concentration (mmol N / m³)",
    )
    push!(axs, ax)
    lines!(ax, times / day, getproperty(timeseries, sym); linewidth=3)
end

fig
```
