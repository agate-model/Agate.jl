# # [Example on the effect of plankton size structures] (@id size_structure_example)

# Install basic modules to run the model
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
nothing #hide

# In this example, we look at how to change the number of phytoplankton and zooplankon size groups given the allometric relationships provided by construct_size_structured_NPZD constructor.

# ## Model with 2 phytoplankton and 2 zooplankton (default NPZD)
# We configure the default plankton model with two phytoplankton and two zooplankton, each represented by different sizes.
phyto_2_zoo_2 = construct_size_structured_NPZD()
nothing #hide

# After generating the plankton model, we can look at the emergent plankton size structure.
model = phyto_2_zoo_2()  # Store the model to avoid multiple function calls
nothing #hide
plankton = ["P1", "P2", "Z1", "Z2"]

# Loop through the plankton and print their diameters
for p in plankton
    println("$(p) diameter: ", model.diameters[p], "um")
end

# To run the model, we need to wrap it into the OceanBioME biogeochemistry model (here light)
bgc_model_phyto_2_zoo_2 = Biogeochemistry(
    phyto_2_zoo_2();
    light_attenuation=FunctionFieldPAR(; grid=BoxModelGrid()), # more intutive if removed (use default) (?)
)
nothing #hide

# Once the plankton model is combined with light availability, we need to wrap it in a physical model.
# Here we are using the OceanBioME Boxmodel().
full_model = BoxModel(; biogeochemistry=bgc_model_phyto_2_zoo_2)

# Next we need to define the initial tracer concentrations.
# "N" states for nutrient (mmol N m-3);
# "D" states for detritus (mmol N m-3);
# "Pi" states for phytoplankton number i (mmol N m-3);
# "Zi" states for zooplankton number i (mmol N m-3).
# Set IC for N and D
params = Dict(:N => 7.0, :D => 0.0)
# Set IC for P1 to P2 and Z1 to Z2 
merge!(params, Dict(Symbol("P$i") => 0.01 for i in 1:2))
merge!(params, Dict(Symbol("Z$i") => 0.05 for i in 1:2))
# Apply parameters to model
set!(full_model; params...)

# ## Simulate model "phyto_2_zoo_2"
# Now that we have setup the model we can run it.
# First, we define the time step (Δt) and the run duration (stop_time).
simulation = Simulation(full_model; Δt=10minutes, stop_time=3years)
nothing #hide

# Next, we define the location where the model outputs are saved:
filename = "box.jld2"
nothing #hide

# In addition, we define the frequence of the model outputs (schedule) and which fields to export:
simulation.output_writers[:fields] = JLD2OutputWriter(
    full_model,
    full_model.fields;
    filename=filename,
    schedule=TimeInterval(1day),
    overwrite_existing=true,
)
# Finally we run the simulation:
run!(simulation)

# ## Plotting phyto_2_zoo_2 model outputs
# Now that we have ran our simulation we can plot it.

# First we need to define the time axis:
timeseries = NamedTuple{keys(full_model.fields)}(
    FieldTimeSeries(filename, "$field")[1, 1, 1, :] for field in keys(full_model.fields)
)

# # Total plankton plots  ## TO FIX, not sure how to add fields together
# # Adding Total Phytoplankton and Total Zooplankton to the timeseries
# timeseries[:TotalP] = timeseries[:P1] .+ timeseries[:P2]  
# timeseries[:TotalZ] = timeseries[:Z1] .+ timeseries[:Z2] 

# # Fields to plot
# fields = [:TotalP, :TotalZ, :D, :N]
# titles = [
#     "Total Phytoplankton"
#     "Total Zooplankton"
#     "Detritus",
#     "Nutrient",
# ]
# nothing #hide

# # Next plot it:
# plots = [
#     plot(; title=titles[i], xlabel="Time (days)", ylabel=string(fields[i])) for
#     i in 1:length(fields)
# ]

# for (i, field) in enumerate(fields)
#     plot!(plots[i], timeseries[field]; color=:blue)
# end

# # Combine
# p = plot(plots...; layout=(2, 2), size=(900, 600))

# savefig(p, "total_phyto_2_zoo_2.png")
# p

# ## Individual plankton plots
# Fields to plot:
fields = [:P1, :P2, :Z1, :Z2]
titles = [
    "Phytoplankton 1",
    "Phytoplankton 2",
    "Zooplankton 1",
    "Zooplankton 2",
]
nothing #hide

# Next plot it:
plots = [
    plot(; title=titles[i], xlabel="Time (days)", ylabel=string(fields[i])) for
    i in 1:length(fields)
]

for (i, field) in enumerate(fields)
    plot!(plots[i], timeseries[field]; color=:blue)
end

p = plot(plots...; layout=(2, 2), size=(900, 600))
p

# # ## Model with 1 phytoplankton and 1 zooplankton (Simple NPZD)

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

# # ## Simulate

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

# # ## Plotting

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

# ## Model with 10 phytoplankton and 10 zooplankton (Complex NPZD)
# Define the model to have 10 phytoplankton and 10 zooplankton
phyto_10_zoo_10 = construct_size_structured_NPZD(; n_phyto=10, n_zoo=10)

# Print emergent parameter values
model = phyto_10_zoo_10()  # Store the model to avoid multiple function calls
plankton = vcat(["P$i" for i in 1:10], ["Z$i" for i in 1:10])
for p in plankton
    println("$(p) diameter: ", model.diameters[p], "um")
end

#Wrap model to the OceanBioME biogeochemistry model:
bgc_model_phyto_10_zoo_10 = Biogeochemistry(
    phyto_10_zoo_10();
    light_attenuation=FunctionFieldPAR(; grid=BoxModelGrid()), # more intutive if removed (use default) (?)
)

# Wrap model to OceanBioME Boxmodel() physical model:
full_model = BoxModel(; biogeochemistry=bgc_model_phyto_10_zoo_10)

# Define initial tracer concentrations
# Set IC for N and D
params = Dict(:N => 7.0, :D => 0.0)
# Set IC for P1 to P10 and Z1 to Z10 
merge!(params, Dict(Symbol("P$i") => 0.01 for i in 1:10))
merge!(params, Dict(Symbol("Z$i") => 0.05 for i in 1:10))
nothing #hide
# Apply parameters to model
set!(full_model; params...)

# ## Run phyto_10_zoo_10 model
# Configure model simulation parameters (time step, duration) and output parameters (save location and frequency):
simulation = Simulation(full_model; Δt=10minutes, stop_time=2years)
filename = "box.jld2"
simulation.output_writers[:fields] = JLD2OutputWriter(
    full_model,
    full_model.fields;
    filename=filename,
    schedule=TimeInterval(1day),
    overwrite_existing=true,
)
run!(simulation)

# ## Plot phyto_10_zoo_10 model outputs

# Define plotting parameters:
timeseries = NamedTuple{keys(full_model.fields)}(
    FieldTimeSeries(filename, "$field")[1, 1, 1, :] for field in keys(full_model.fields)
)

# Plot phytoplankton
fields = vcat([Symbol("P$i") for i in 1:10])
titles = vcat(
    ["Phytoplankton $i" for i in 1:10],  # Phytoplankton 1 to 10
)
plots = [
    plot(; title=titles[i], xlabel="Time (days)", ylabel=string(fields[i])) for
    i in 1:length(fields)
]
for (i, field) in enumerate(fields)
    plot!(plots[i], timeseries[field]; color=:blue)
end
p = plot(plots...; layout=(5, 2), size=(900, 600))
p

# Plot zooplankton
fields = vcat([Symbol("Z$i") for i in 1:10])
titles = vcat(
    ["Zooplankton $i" for i in 1:10],  # Zooplankton 1 to 10
)
plots = [
    plot(; title=titles[i], xlabel="Time (days)", ylabel=string(fields[i])) for
    i in 1:length(fields)
]
for (i, field) in enumerate(fields)
    plot!(plots[i], timeseries[field]; color=:blue)
end
p = plot(plots...; layout=(5, 2), size=(900, 600))
p

# Plot others
fields = [:D, :N]
titles = [
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
p = plot(plots...; layout=(1, 2), size=(900, 600))
p