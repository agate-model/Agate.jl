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

# In this example, we look at how to increase the number of phytoplankton and zooplankton size groups based on the allometric relationships defined in the construct_size_structured_NPZD constructor. 
# Allometric relationships describe how biological rates (e.g., growth, grazing, and mortality) scale with organism size. These relationships serve as fundamental principles in trait-based ecosystem modeling, where individual traits like size influence ecosystem-scale dynamics, allowing models to capture the diversity and functionality of plankton communities. 
# Size is a particularly critical trait, as it strongly influences individual physiology, ecological interactions, and environmental adaptations. Smaller phytoplankton, for example, have a higher surface-area-to-volume ratio, which enhances their nutrient uptake efficiency. This makes them well-suited for survival in nutrient-poor environments, where resource acquisition is a key limiting factor. 

# Here, we give the example of modelling an ecosystem with 4 phytoplankton and 4 zooplankton size classes.

# Allometric relationships usually scale as a power law: M = aV^b 


# ## Contruct the model
# Define the model to have 4 phytoplankton and 4 zooplankton size classes with a range of 6-150 \mum and 60-1500 \mum, respectively, distributed equally on a log scare.
phyto_4_zoo_4 = construct_size_structured_NPZD(; 
n_phyto=4, 
n_zoo=4,
phyto_diameters=Dict(
"min_diameter" => 6, "max_diameter" => 150, "splitting" => "log_splitting"
),
zoo_diameters=Dict(
    "min_diameter" => 60, "max_diameter" => 1500, "splitting" => "log_splitting"
))
nothing #hide

# Print emergent parameter values
model = phyto_4_zoo_4()  # Store the model to avoid multiple function calls
plankton = vcat(["P$i" for i in 1:4], ["Z$i" for i in 1:4])
for p in plankton
    println("$(p) diameter: ", model.diameters[p], "um")
end

# Wrap model to the OceanBioME biogeochemistry model:
bgc_model_phyto_4_zoo_4 = Biogeochemistry(
    phyto_4_zoo_4();
    light_attenuation=FunctionFieldPAR(; grid=BoxModelGrid()), # more intutive if removed (use default) (?)
)
nothing #hide

# Wrap model to OceanBioME Boxmodel() physical model:
full_model = BoxModel(; biogeochemistry=bgc_model_phyto_4_zoo_4)
nothing #hide

# Define initial tracer concentrations
# Set Initial Conditions for N and D
params = Dict(:N => 7.0, :D => 0.0)
# Set IC for P4 to P5 and Z1 to Z4 
merge!(params, Dict(Symbol("P$i") => 0.01 for i in 1:4))
merge!(params, Dict(Symbol("Z$i") => 0.05 for i in 1:4))
nothing #hide
# Apply parameters to model
set!(full_model; params...)

# ## Run model
# Configure model simulation parameters (time step, duration) and output parameters (save location and frequency):
simulation = Simulation(full_model; Δt=10minutes, stop_time=4years)
filename = "box.jld2"
simulation.output_writers[:fields] = JLD2OutputWriter(
    full_model,
    full_model.fields;
    filename=filename,
    schedule=TimeInterval(1day),
    overwrite_existing=true,
)
run!(simulation)
nothing #hide

# ## Plot model outputs

# Define plotting parameters:
timeseries = NamedTuple{keys(full_model.fields)}(
    FieldTimeSeries(filename, "$field")[1, 1, 1, :] for field in keys(full_model.fields)
)
nothing #hide

# Plot phytoplankton
fields = vcat([Symbol("P$i") for i in 1:4])
titles = vcat(
    ["Phytoplankton $i" for i in 1:4],  # Phytoplankton 1 to 4
)
plots = [
    plot(; title=titles[i], xlabel="Time (days)", ylabel=string(fields[i])) for
    i in 1:length(fields)
]
for (i, field) in enumerate(fields)
    plot!(plots[i], timeseries[field]; color=:blue)
end
p1 = plot(plots...; layout=(2, 2), size=(800, 600))
p1

# Plot zooplankton
fields = vcat([Symbol("Z$i") for i in 1:4])
titles = vcat(
    ["Zooplankton $i" for i in 1:4],  # Zooplankton 1 to 4
)
plots = [
    plot(; title=titles[i], xlabel="Time (days)", ylabel=string(fields[i])) for
    i in 1:length(fields)
]
for (i, field) in enumerate(fields)
    plot!(plots[i], timeseries[field]; color=:blue)
end
p2 = plot(plots...; layout=(2, 2), size=(800, 600))
p2

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
p3 = plot(plots...; layout=(1, 2), size=(900, 600))
p3