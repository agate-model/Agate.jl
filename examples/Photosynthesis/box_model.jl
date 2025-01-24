using Agate
using Agate.Library.Light
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans.Units
using Plots
using Agate.Models.Constructors

const year = years = 365day

# ==================================================
# Define BGC model 
# ==================================================

N2P2ZD_default_photosynthesis = construct_size_structured_NPZD()

bgc_model_default_photosynthesis = Biogeochemistry(
    N2P2ZD_default_photosynthesis();
    light_attenuation=FunctionFieldPAR(; grid=BoxModelGrid()),
)

full_model_default_photosynthesis = BoxModel(;
    biogeochemistry=bgc_model_default_photosynthesis
)
set!(full_model_default_photosynthesis; P1=0.01, P2=0.01, Z1=0.05, Z2=0.05, N=7.0, D=1)

# ==================================================
# Simulate
# ==================================================

filename_default_photosynthesis = "box_n2p2zd_default_photosynthesis.jld2"

simulation_default_photosynthesis = Simulation(
    full_model_default_photosynthesis; Δt=10minutes, stop_time=3years
)
simulation_default_photosynthesis.output_writers[:fields] = JLD2OutputWriter(
    full_model_default_photosynthesis,
    full_model_default_photosynthesis.fields;
    filename=filename,
    schedule=TimeInterval(1day),
    overwrite_existing=true,
)

run!(simulation_default_photosynthesis)

timeseries_default_photosynthesis = NamedTuple{
    keys(full_model_default_photosynthesis.fields)
}(
    FieldTimeSeries(filename_default_photosynthesis, "$field")[1, 1, 1, :] for
    field in keys(full_model_default_photosynthesis.fields)
)

# ==================================================
# Plotting
# ==================================================

p1 = plot(timeseries_default_photosynthesis.P1; label="P1", title="Phytoplankton 1")
p2 = plot(timeseries_default_photosynthesis.P2; label="P2", title="Phytoplankton 2")
p3 = plot(timeseries_default_photosynthesis.Z1; label="Z1", title="Zooplankton 1")
p4 = plot(timeseries_default_photosynthesis.Z2; label="Z2", title="Zooplankton 2")
p5 = plot(timeseries_default_photosynthesis.D; label="D", title="Detritus")
p6 = plot(timeseries_default_photosynthesis.N; label="N", title="nutrient")

# Arrange plots in a 2x3 layout
p = plot(p1, p2, p3, p4, p5, p6; layout=(2, 3), size=(900, 600))

savefig(p, "N2P2ZD_box_default_photosynthesis.png")
