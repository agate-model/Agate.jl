# This example shows how to create a Biogeochemistry (BGC) model from parameter and tracer
# definitions and run a box model (0D) simulation using OceanBioME and Oceananigans.

using Agate
using Agate.Library.Light
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans.Units
using Plots

const year = years = 365day

# ==================================================
# Define BGC model (NPZD with cyclical PAR)
# ==================================================

include(joinpath("NPZD", "tracers.jl"))
bgc_model = Biogeochemistry(
    NPZD(); light_attenuation=FunctionFieldPAR(; grid=BoxModelGrid())
)
full_model = BoxModel(; biogeochemistry=bgc_model)
set!(full_model; N=7.0, P=0.01, Z=0.05, D=0.0)

# ==================================================
# Simulate
# ==================================================

filename = "box.jld2"

simulation = Simulation(full_model; Δt=5minutes, stop_time=3years)
simulation.output_writers[:fields] = JLD2OutputWriter(
    full_model,
    full_model.fields;
    filename=filename,
    schedule=TimeInterval(1day),
    overwrite_existing=true,
)

run!(simulation)

timeseries = NamedTuple{keys(full_model.fields)}(
    FieldTimeSeries(filename, "$field")[1, 1, 1, :] for field in keys(full_model.fields)
)

# ==================================================
# Plotting
# ==================================================

p = plot(timeseries.P; label="P")
plot!(p, timeseries.Z; label="Z")
plot!(p, timeseries.D; label="D")
savefig(p, "NPZD_box_oceananigans.png")
