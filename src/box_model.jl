using OceanBioME
using Oceananigans.Fields: FunctionField
using Oceananigans.Units

export run_npzd_boxmodel

const year = years = 365day

# do we ever want to change this constant?
const z = -10 # specify the nominal depth of the box for the PAR profile

# should this be defined somewhere else in the library?
PAR⁰(t) = 60 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2
PAR_f(t) = PAR⁰(t) * exp(0.2z) # Modify the PAR based on the nominal depth and exponential decay


function run_npzd_boxmodel(init_cond, parameters; Δt=5minutes, stop_time=3years, save_interval=8hours, filename="box.jld2", overwrite=true)

    N,P,Z,D = init_cond
    # are there other parameters we want to pass here?
    base_maximum_growth, maximum_grazing_rate = parameters

    grid = BoxModelGrid() # 1x1x1 grid
    clock = Clock(time = zero(grid))
    PAR = FunctionField{Center, Center, Center}(PAR_f, grid; clock)

    biogeochemistry = NutrientPhytoplanktonZooplanktonDetritus(;
        grid,
        base_maximum_growth = base_maximum_growth,
        maximum_grazing_rate = maximum_grazing_rate,
        # this just returns the value of PAR at time t
        light_attenuation_model = PrescribedPhotosyntheticallyActiveRadiation(PAR),
        # by default P and D are set to sink at a constant rate although it probably
        # doesn't matter whether this is set or not for a BoxModel (unless it has an
        # open bottom?)
        sinking_speeds = NamedTuple()
    )

    # temperature is set to 0 by default
    model = BoxModel(;biogeochemistry, clock)
    # to access value of each tracer: model.fields.$tracer.data[1,1,1]
    set!(model, N=N, P=P, Z=Z, D=D)

    simulation = Simulation(model; Δt = Δt, stop_time = stop_time)
    simulation.output_writers[:fields] = JLD2OutputWriter(model,
                            model.fields;
                            filename = filename,
                            schedule = TimeInterval(save_interval),
                            overwrite_existing = overwrite)
    run!(simulation)

    timeseries = NamedTuple{keys(model.fields)}(FieldTimeSeries("box.jld2", "$field")[1, 1, 1, :] for field in keys(model.fields))

    return timeseries
end
