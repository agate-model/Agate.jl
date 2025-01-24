using Agate
using Agate.Library.Light
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans.Units
using Plots
using Agate.Models.Constructors
using Agate.Models.Tracers
using Agate.Library.Photosynthesis

const year = years = 365day

# ==================================================
# Define BGC model 
# ==================================================

N2P2ZD_geider_photosynthesis = construct_size_structured_NPZD(;
    phyto_args=Dict(
        "diameters" =>
            Dict("min_diameter" => 2, "max_diameter" => 10, "splitting" => "log_splitting"),
        "allometry" => Dict(
            "maximum_growth_rate" => Dict("a" => 2 / day, "b" => -0.15),
            "nutrient_half_saturation" => Dict("a" => 0.17, "b" => 0.27),
        ),
        "linear_mortality" => 8e-7 / second,
        "photosynthetic_slope" => 0.3953 / day,
        "chlorophyll_to_carbon_ratio" => 3,
    ),
    nutrient_dynamics=nutrients_geider_light,
    phyto_dynamics=phytoplankton_growth_single_nutrient_geider_light,
)

bgc_model_geider_photosynthesis = Biogeochemistry(
    N2P2ZD_geider_photosynthesis();
    light_attenuation=FunctionFieldPAR(; grid=BoxModelGrid()),
)

full_model_geider_photosynthesis = BoxModel(;
    biogeochemistry=bgc_model_geider_photosynthesis
)
set!(full_model_geider_photosynthesis; P1=0.01, P2=0.01, Z1=0.05, Z2=0.05, N=7.0, D=1)

# ==================================================
# Simulate
# ==================================================

filename_geider_photosynthesis = "box_n2p2zd_geider_photosynthesis.jld2"

simulation_geider_photosynthesis = Simulation(
    full_model_geider_photosynthesis; Δt=10minutes, stop_time=3years
)
simulation_geider_photosynthesis.output_writers[:fields] = JLD2OutputWriter(
    full_model_geider_photosynthesis,
    full_model_geider_photosynthesis.fields;
    filename=filename_geider_photosynthesis,
    schedule=TimeInterval(1day),
    overwrite_existing=true,
)

run!(simulation_geider_photosynthesis)

timeseries_geider_photosynthesis = NamedTuple{keys(full_model_geider_photosynthesis.fields)}(
    FieldTimeSeries(filename_geider_photosynthesis, "$field")[1, 1, 1, :] for
    field in keys(full_model_geider_photosynthesis.fields)
)

# ==================================================
# Plotting
# ==================================================

p1 = plot(timeseries_geider_photosynthesis.P1; label="P1", title="Phytoplankton 1")
p2 = plot(timeseries_geider_photosynthesis.P2; label="P2", title="Phytoplankton 2")
p3 = plot(timeseries_geider_photosynthesis.Z1; label="Z1", title="Zooplankton 1")
p4 = plot(timeseries_geider_photosynthesis.Z2; label="Z2", title="Zooplankton 2")
p5 = plot(timeseries_geider_photosynthesis.D; label="D", title="Detritus")
p6 = plot(timeseries_geider_photosynthesis.N; label="N", title="nutrient")

# Arrange plots in a 2x3 layout
p = plot(p1, p2, p3, p4, p5, p6; layout=(2, 3), size=(900, 600))

savefig(p, "N2P2ZD_box_geider_photosynthesis.png")
