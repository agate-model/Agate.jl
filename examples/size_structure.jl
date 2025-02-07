using Agate
using Agate.Library.Light
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans.Units
using Plots
using Agate.Constructors
using Agate.Models.Tracers
using Agate.Library.Photosynthesis

const year = years = 365day

# Some description

# the default model will define two phytoplankton and two zooplankton
phyto_2_zoo_2 = construct_size_structured_NPZD()

# After generating the model we can look at the emergent parameter values
println("Phytoplankton 1 diameter: ", phyto_2_zoo_2().diameters["P1"], "um")
println("Phytoplankton 2 diameter: ", phyto_2_zoo_2().diameters["P2"], "um")
println("Zooplankton 1 diameter: ", phyto_2_zoo_2().diameters["Z1"], "um")
println("Zooplankton 2 diameter: ", phyto_2_zoo_2().diameters["Z2"], "um")

#to change the number of phytoplankton we can specify n_phyto or n_zoo
phyto_1_zoo_1 = construct_size_structured_NPZD(; n_phyto=1, n_zoo=1)

#to run the model we need to wrap it into a OceanBioME biogeochemistry model
bgc_model_phyto_2_zoo_2 = Biogeochemistry(
    phyto_2_zoo_2();
    light_attenuation=FunctionFieldPAR(; grid=BoxModelGrid()), # more intutive if removed (use default) (?)
)

# once we have our combined light and ecosystem model we need to define the physical setting
# here we are using a OceanBioME Boxmodel()
full_model = BoxModel(; biogeochemistry=bgc_model_phyto_2_zoo_2)

# Next we need to define the initial tracer concentrations.
# Here "N" is a nutrient (mmol N m-3)
# "D" is detritus (mmol N m-3)
# "Pi" are the phytoplankton (mmol N m-3)
# "Zi" are the zooplankton (mmol N m-3)
set!(full_model; N=7.0, P1=0.01, P2=0.01, Z1=0.05, Z2=0.05, D=0.0)

# Now that we have setup the model we can run it
# First we need to define the save location:
filename = "box.jld2"
# Next we define the time step (Δt) and the run duration (stop_time)
simulation = Simulation(full_model; Δt=5minutes, stop_time=3years)
# In addition we need to define how frequently to save the model (schedule=TimeInterval(1day)) and which fields to export
simulation.output_writers[:fields] = JLD2OutputWriter(
    full_model,
    full_model.fields;
    filename=filename,
    schedule=TimeInterval(1day),
    overwrite_existing=true,
)
# Finally we run the simulation:
run!(simulation)

# Now that we have ran our simulation we can plot it:

# First we need to define the time axis:
timeseries = NamedTuple{keys(full_model.fields)}(
    FieldTimeSeries(filename, "$field")[1, 1, 1, :] for field in keys(full_model.fields)
)

# Next plot it:

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
    plot!(plots[i], timeseries[field]; color=:blue)
end

p = plot(plots...; layout=(2, 3), size=(900, 600))

savefig(p, "2_phyto_2_zoo.png")
