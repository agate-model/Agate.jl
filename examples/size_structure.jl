# Read the basic modules to run the model
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

# In this example, we look at how to change the number of phytoplankton and zooplankon size groups, given the allometric relationships provided by construct_size_structured_NPZD constructor

# ==================================================
# Define model with 2 phytoplankton and 2 zooplankton (default NPZD)
# ==================================================

# We define the default model to have two phytoplankton and two zooplankton (characterised with different size)
phyto_2_zoo_2 = construct_size_structured_NPZD()

# After generating the model, we can look at the emergent parameter values
model = phyto_2_zoo_2()  # Store the model to avoid multiple function calls

# Define the keys in an ordered array
plankton = ["P1", "P2", "Z1", "Z2"]

# Loop through the plankton and print their diameters
for p in plankton
    println("$(p) diameter: ", model.diameters[p], "um")
end
# println("Phytoplankton 1 diameter: ", phyto_2_zoo_2().diameters["P1"], "um")
# println("Phytoplankton 2 diameter: ", phyto_2_zoo_2().diameters["P2"], "um")
# println("Zooplankton 1 diameter: ", phyto_2_zoo_2().diameters["Z1"], "um")
# println("Zooplankton 2 diameter: ", phyto_2_zoo_2().diameters["Z2"], "um")

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


# ==================================================
# Simulate
# ==================================================
# Now that we have setup the model we can run it
# First, we define the time step (Δt) and the run duration (stop_time)
simulation = Simulation(full_model; Δt=5minutes, stop_time=3years)
# Next, we define the location where to save the model outputs:
filename = "box.jld2"
# In addition, we define the frequence of the model outputs (schedule) and which fields to export
simulation.output_writers[:fields] = JLD2OutputWriter(
    full_model,
    full_model.fields;
    filename=filename,
    schedule=TimeInterval(1day),
    overwrite_existing=true,
)
# Finally we run the simulation:
run!(simulation)

# ==================================================
# Plotting
# ==================================================
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

savefig(p, "phyto_2_zoo_2.png")

# # ==================================================
# # Define model with 1 phytoplankton and 1 zooplankton (Simple NPZD)
# # ==================================================

# # Define the model to have one phytoplankton and one zooplankton
# phyto_1_zoo_1 = construct_size_structured_NPZD(; n_phyto=1, n_zoo=1)

# # Print emergent parameter values
# model = phyto_1_zoo_1()  # Store the model to avoid multiple function calls
# plankton = ["P1", "Z1"]
# for p in plankton
#     println("$(p) diameter: ", model.diameters[p], "um")
# end

# #Wrap model to the OceanBioME biogeochemistry model
# bgc_model_phyto_1_zoo_1 = Biogeochemistry(
#     phyto_1_zoo_1();
#     light_attenuation=FunctionFieldPAR(; grid=BoxModelGrid()), # more intutive if removed (use default) (?)
# )

# # Wrap model to the physical model
# # here we are using a OceanBioME Boxmodel()
# full_model = BoxModel(; biogeochemistry=bgc_model_phyto_1_zoo_1)

# # Define initial tracer concentrations
# # Here "N" is a nutrient (mmol N m-3)
# # "D" is detritus (mmol N m-3)
# # "Pi" are the phytoplankton (mmol N m-3)
# # "Zi" are the zooplankton (mmol N m-3)
# set!(full_model; N=7.0, P1=0.01, Z1=0.05, D=0.0)

# # ==================================================
# # Simulate
# # ==================================================

# # Configure model simulation parameters (time step, duration) and output parameters (save location and frequency):
# simulation = Simulation(full_model; Δt=5minutes, stop_time=3years)
# filename = "box.jld2"
# simulation.output_writers[:fields] = JLD2OutputWriter(
#     full_model,
#     full_model.fields;
#     filename=filename,
#     schedule=TimeInterval(1day),
#     overwrite_existing=true,
# )

# run!(simulation)

# # ==================================================
# # Plotting
# # ==================================================

# # Define plotting parameters:
# timeseries = NamedTuple{keys(full_model.fields)}(
#     FieldTimeSeries(filename, "$field")[1, 1, 1, :] for field in keys(full_model.fields)
# )
# fields = [:P1, :Z1, :D, :N]
# titles = [
#     "Phytoplankton 1",
#     "Zooplankton 1",
#     "Detritus",
#     "Nutrient",
# ]

# plots = [
#     plot(; title=titles[i], xlabel="Time (days)", ylabel=string(fields[i])) for
#     i in 1:length(fields)
# ]

# for (i, field) in enumerate(fields)
#     plot!(plots[i], timeseries[field]; color=:blue)
# end

# p = plot(plots...; layout=(2, 2), size=(900, 600))

# savefig(p, "phyto_1_zoo_1.png")

# ==================================================
# Define model with 10 phytoplankton and 10 zooplankton (Complex NPZD)
# ==================================================

# Define the model to have 10 phytoplankton and 10 zooplankton
phyto_10_zoo_10 = construct_size_structured_NPZD(; n_phyto=10, n_zoo=10)

# Print emergent parameter values
model = phyto_10_zoo_10()  # Store the model to avoid multiple function calls
plankton = vcat(["P$i" for i in 1:10], ["Z$i" for i in 1:10])
for p in plankton
    println("$(p) diameter: ", model.diameters[p], "um")
end

#Wrap model to the OceanBioME biogeochemistry model
bgc_model_phyto_10_zoo_10 = Biogeochemistry(
    phyto_10_zoo_10();
    light_attenuation=FunctionFieldPAR(; grid=BoxModelGrid()), # more intutive if removed (use default) (?)
)

# Wrap model to the physical model
# here we are using a OceanBioME Boxmodel()
full_model = BoxModel(; biogeochemistry=bgc_model_phyto_10_zoo_10)

# Define initial tracer concentrations
# Here "N" is a nutrient (mmol N m-3)
# "D" is detritus (mmol N m-3)
# "Pi" are the phytoplankton (mmol N m-3)
# "Zi" are the zooplankton (mmol N m-3)
set!(full_model; N=7.0, P1=0.01, Z1=0.05, D=0.0)

# ==================================================
# Simulate
# ==================================================

# Configure model simulation parameters (time step, duration) and output parameters (save location and frequency):
simulation = Simulation(full_model; Δt=5minutes, stop_time=3years)
filename = "box.jld2"
simulation.output_writers[:fields] = JLD2OutputWriter(
    full_model,
    full_model.fields;
    filename=filename,
    schedule=TimeInterval(1day),
    overwrite_existing=true,
)

run!(simulation)

# ==================================================
# Plotting
# ==================================================

# Define plotting parameters:
timeseries = NamedTuple{keys(full_model.fields)}(
    FieldTimeSeries(filename, "$field")[1, 1, 1, :] for field in keys(full_model.fields)
)
fields = [:P1, :Z1, :D, :N]
titles = [
    "Phytoplankton 1",
    "Zooplankton 1",
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

p = plot(plots...; layout=(2, 2), size=(900, 600))

savefig(p, "phyto_10_zoo_10.png")
