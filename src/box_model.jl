include("Library/light.jl")

using .Light

using OceanBioME
using Oceananigans
using Oceananigans: Clock

using Oceananigans.Units
using Oceananigans.Fields: FunctionField

import Oceananigans: set!

export create_box_model, run_box_model

const year = years = 365day

"""
    create_box_model(
        bgc_model, init_conditions; PAR_f=cyclical_PAR, PAR_parameters=(; z=-10)
    ) -> OceanBioME.BoxModel

Create an OceanBioME.BoxModel object and set initial values.

# Arguments
- `bgc_model`: biogeochemistry model, a subtype of AbstractContinuousFormBiogeochemistry,
    e.g., returned by `Agate.create_bgc_struct()`
- `init_conditions`: NamedTuple of initial values

# Keywords
- `PAR_f`: a time dependant PAR function (defaults to `Agate.Library.Light.cyclical_PAR`)
- `PAR_parameters`: any fixed parameters of the PAR function (e.g., depth `z`) passed as NamedTuple
"""
function create_box_model(
    bgc_model, init_conditions; PAR_f=cyclical_PAR, PAR_parameters=(; z=-10)
)
    grid = BoxModelGrid() # 1x1x1 grid
    clock = Clock(; time=zero(grid))
    if isnothing(PAR_f)
        light_attenuation = nothing
    else
        PAR = FunctionField{Center,Center,Center}(
            PAR_f, grid; clock, parameters=PAR_parameters
        )
        light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(PAR)
    end

    biogeochemistry = Biogeochemistry(bgc_model; light_attenuation=light_attenuation)

    model = BoxModel(; biogeochemistry, clock)
    set!(model, init_conditions)

    return model
end

"""
    run_box_model(bgc_model, init_conditions; kwargs...) -> NamedTuple

Returns timeseries for each tracer of the form (<tracer name>: [<value at t1>, ...], ...)
(results are also saved to a file).

# Arguments
- `bgc_model`: biogeochemistry model, a subtype of AbstractContinuousFormBiogeochemistry
    (e.g., returned by `Agate.create_bgc_struct`)
- `init_conditions`: NamedTuple of initial values

# Keywords
- `PAR_f`: a time dependant PAR function (defaults to `Agate.Library.Light.cyclical_PAR`)
- `PAR_parameters`: any fixed parameters of the PAR function (e.g., depth `z`) passed as NamedTuple
- `Δt``: simulation step time
- `stop_time`: until when to run the simulation
- `save_interval`: interval at which to save simulation results
- `filename`: name of file to save simulation results to
- `overwrite`: whether to overwrite existing files
"""
function run_box_model(
    bgc_model,
    init_conditions;
    PAR_f=cyclical_PAR,
    PAR_parameters=(; z=-10),
    Δt=5minutes,
    stop_time=3years,
    save_interval=1day,
    filename="box.jld2",
    overwrite=true,
)
    model = create_box_model(
        bgc_model, init_conditions; PAR_f=PAR_f, PAR_parameters=PAR_parameters
    )

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
        FieldTimeSeries("box.jld2", "$field")[1, 1, 1, :] for field in keys(model.fields)
    )

    return timeseries
end

"""
    set!(model::BoxModel, init_conditions::NamedTuple)

Set the `BoxModel` initial conditions (Field values).

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
