using Agate
using Agate.Library.Light
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans.Units
using Plots

const year = years = 365day

include(joinpath("N2P2ZD", "tracers.jl"))

using Adapt, CUDA
CUDA.allowscalar(true)

# Adapt.@adapt_structure N2P2ZD


adapted_instance = Adapt.adapt(CuArray, N2P2ZD())


bgc_model = Biogeochemistry(
    adapted_instance; light_attenuation=FunctionFieldPAR(; grid=BoxModelGrid(arch=GPU()))
)
full_model = BoxModel(; biogeochemistry=bgc_model)
set!(full_model; N=7.0, P1=0.01, Z1=0.05, P2=0.01, Z2=0.05, D=0.0)

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

p = plot(timeseries.P1; label="P1")
plot!(p, timeseries.Z1; label="Z1")
plot!(p, timeseries.D; label="D")
savefig(p, "N2P2ZD_box_GPU_manual_adapt.png")
