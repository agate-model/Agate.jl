using Agate
using Agate.Library.Light
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans.Units
using Plots
using Agate.Constructors: NPZD_size_structured
using Agate.Models.Tracers
using Agate.Library.Photosynthesis

const year = years = 365day

# ==================================================
# Define BGC models
# ==================================================

# Default photosynthesis model
N2P2ZD_default_photosynthesis = NPZD_size_structured.construct()
bgc_model_default_photosynthesis = Biogeochemistry(
    N2P2ZD_default_photosynthesis();
    light_attenuation=FunctionFieldPAR(; grid=BoxModelGrid()),
)
full_model_default_photosynthesis = BoxModel(;
    biogeochemistry=bgc_model_default_photosynthesis
)
set!(full_model_default_photosynthesis; P1=0.01, P2=0.01, Z1=0.05, Z2=0.05, N=7.0, D=1)

# Geider photosynthesis model
N2P2ZD_geider = NPZD_size_structured.construct(;
    phyto_args=NPZD_size_structured.DEFAULT_PHYTO_GEIDER_ARGS,
    nutrient_dynamics=nutrients_geider_light,
    phyto_dynamics=phytoplankton_growth_single_nutrient_geider_light,
)
N2P2ZD_geider_photosynthesis = NPZD_size_structured.instantiate(
    N2P2ZD_geider; phyto_args=NPZD_size_structured.DEFAULT_PHYTO_GEIDER_ARGS
)
bgc_model_geider_photosynthesis = Biogeochemistry(
    N2P2ZD_geider_photosynthesis; light_attenuation=FunctionFieldPAR(; grid=BoxModelGrid())
)
full_model_geider_photosynthesis = BoxModel(;
    biogeochemistry=bgc_model_geider_photosynthesis
)
set!(full_model_geider_photosynthesis; P1=0.01, P2=0.01, Z1=0.05, Z2=0.05, N=7.0, D=1)

# ==================================================
# Simulate models
# ==================================================

function run_simulation!(model, filename)
    simulation = Simulation(model; Δt=10minutes, stop_time=3years)
    simulation.output_writers[:fields] = JLD2OutputWriter(
        model,
        model.fields;
        filename=filename,
        schedule=TimeInterval(1day),
        overwrite_existing=true,
    )
    run!(simulation)
    return NamedTuple{keys(model.fields)}(
        FieldTimeSeries(filename, "$field")[1, 1, 1, :] for field in keys(model.fields)
    )
end

filename_default_photosynthesis = "box_n2p2zd_default_photosynthesis.jld2"
filename_geider_photosynthesis = "box_n2p2zd_geider_photosynthesis.jld2"

timeseries_default_photosynthesis = run_simulation!(
    full_model_default_photosynthesis, filename_default_photosynthesis
)
timeseries_geider_photosynthesis = run_simulation!(
    full_model_geider_photosynthesis, filename_geider_photosynthesis
)

# ==================================================
# Plotting
# ==================================================

fields = [:P1, :P2, :Z1, :Z2, :D, :N]
titles = [
    "Phytoplankton 1",
    "Phytoplankton 2",
    "Zooplankton 1",
    "Zooplankton 2",
    "Detritus",
    "Nutrient",
]

plots = [
    plot(; title=titles[i], xlabel="Time (days)", ylabel=string(fields[i])) for
    i in 1:length(fields)
]

for (i, field) in enumerate(fields)
    plot!(plots[i], timeseries_default_photosynthesis[field]; label="Default", color=:blue)
    plot!(plots[i], timeseries_geider_photosynthesis[field]; label="Geider", color=:green)
end

p = plot(plots...; layout=(2, 3), size=(900, 600))
savefig(p, "N2P2ZD_contrast_photosynthesis.png")
