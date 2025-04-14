# In this model, a single plankton functional type is defined, which can then be specified as predator or prey based on the parameter values.
# For example, for prey # maximum_grazing_rate are set to 0, while for predators maximum_growth_rate are set to 0.

using Agate
using Agate.Library.Light
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans.Units
using Plots

const year = years = 365day

# ==================================================
# Define BGC model (N2P2ZD)
# ==================================================

include("tracers.jl")
bgc_model = Biogeochemistry(
    N2P2ZD(); light_attenuation=FunctionFieldPAR(; grid=BoxModelGrid())
)
full_model = BoxModel(; biogeochemistry=bgc_model)
set!(full_model; N=7.0, Z2=0.05, D=0.0, P1=0.1, P2=0.1, Z1=0.05)

# ==================================================
# Simulate
# ==================================================

filename = "box_n2p2zd.jld2"

simulation = Simulation(full_model; Δt=10minutes, stop_time=10years)
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

p1 = plot(timeseries.P1; label="P1", title="Phytoplankton 1")
p2 = plot(timeseries.P2; label="P2", title="Phytoplankton 2")
p3 = plot(timeseries.Z1; label="Z1", title="Zooplankton 1")
p4 = plot(timeseries.Z2; label="Z2", title="Zooplankton 2")
p5 = plot(timeseries.D; label="D", title="Detritus")
p6 = plot(timeseries.N; label="N", title="nutrient")

# Arrange plots in a 2x3 layout
p = plot(p1, p2, p3, p4, p5, p6; layout=(2, 3), size=(900, 600))

savefig(p, "N2P2ZD_box.png")
