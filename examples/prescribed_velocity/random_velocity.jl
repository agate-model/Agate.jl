using Oceananigans
using  Oceananigans.Units
using CairoMakie

const year = years = 365day

#set_theme!(Theme(fontsize = 20, linewidth=3))

grid = RectilinearGrid(size=48, z=(-0.5, 0.5), topology=(Flat, Flat, Bounded))


#prescribed physics:

@inline w(z, t) =  randn()

velocities = PrescribedVelocityFields(; w)


model = HydrostaticFreeSurfaceModel(; grid, velocities=velocities, tracers=:T)

width = 0.1
initial_temperature(z) = exp(-z^2 / (2width^2))
set!(model, T=initial_temperature)
# Time-scale for diffusion across a grid cell

simulation = Simulation(model, Δt = 10minutes, stop_iteration = 1day)

save_interval = 10minutes
simulation.output_writers[:temperature] =
    JLD2OutputWriter(model, model.tracers,
                     filename = "one_dimensional_diffusion.jld2",
                     schedule=TimeInterval(save_interval),
                     overwrite_existing = true)

run!(simulation)


T_timeseries = FieldTimeSeries("one_dimensional_diffusion.jld2", "T")
times = T_timeseries.times


data = T_timeseries.data[1, 1, :, :]  # This will give you a depth×time matrix

data = Matrix(data)

print(data)

# Define the axes: depth and time
depth = 1:12
time = 1:101

# Create the heatmap
heatmap(time, depth, data)
