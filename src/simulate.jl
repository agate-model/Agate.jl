using OceanBioME
using Oceananigans

using Oceananigans.Units

import Oceananigans: set!

export run_simulation

const year = years = 365day

"""
    run_simulation(model, init_conditions; kwargs...) -> NamedTuple

Returns timeseries for each tracer of the form (<tracer name>: [<value at t1>, ...], ...)
(results are also saved to a file).

# Arguments
- `model`: BoxModel
- `init_conditions`: NamedTuple of initial values

# Keywords
- `Δt`: simulation step time
- `stop_time`: until when to run the simulation
- `save_interval`: interval at which to save simulation results
- `filename`: name of file to save simulation results to
- `overwrite`: whether to overwrite existing files
"""
function run_simulation(
    model::BoxModel,
    init_conditions;
    Δt=5minutes,
    stop_time=3years,
    save_interval=1day,
    filename="box.jld2",
    overwrite=true,
)
    set!(model, init_conditions)

    simulation = Simulation(model; Δt=Δt, stop_time=stop_time)
    simulation.output_writers[:fields] = JLD2OutputWriter(
        model,
        model.fields;
        filename=filename,
        schedule=TimeInterval(save_interval),
        overwrite_existing=overwrite,
    )
    run!(simulation)

    # to access value of each tracer: model.fields.$tracer.data[1,1,1]
    timeseries = NamedTuple{keys(model.fields)}(
        FieldTimeSeries(filename, "$field")[1, 1, 1, :] for field in keys(model.fields)
    )

    return timeseries
end

"""
    set!(model::BoxModel, init_conditions::NamedTuple)

Set the `BoxModel` initial conditions (tracer values).

# Arguments
- `model`: the model to set the arguments for
- `init_conditions`: NamedTuple of initial values
"""
function set!(model::BoxModel, init_conditions::NamedTuple)
    for (fldname, value) in pairs(init_conditions)
        if fldname ∈ propertynames(model.fields)
            ϕ = getproperty(model.fields, fldname)
        else
            throw(ArgumentError("name $fldname not found in model.fields."))
        end
        set!(ϕ, value)
    end
    return nothing
end
