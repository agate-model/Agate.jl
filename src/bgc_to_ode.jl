using OrdinaryDiffEq

using Oceananigans.Biogeochemistry:
    required_biogeochemical_tracers, required_biogeochemical_auxiliary_fields

"""
    bgc_to_ode(biogeochemistry, PAR_f, init_conditions, tspan, parameters=nothing) -> ODEProblem

Generate ODEProblem from Biogeochemistry (suitable for box models).

# Arguments
- `biogeochemistry`: subtype of Oceananigans.Biogeochemistry
- `PAR_f`: a time dependant PAR function
- `init_conditions`: a NamedTuple of initial tracer values
- `tspan`: a Tuple of Floats indicating simulation start and stop times
- `parameters`: ordered names of parameters that will be pasesed to ODEProblem, specified as
   a Tuple of Symbols. Defaults to `nothing`, in which case all `biogeochemistry` parameters
   are used.
"""
function bgc_to_ode(biogeochemistry, PAR_f, init_conditions, tspan, parameters=nothing)

    # retrieve tracer order
    model = biogeochemistry()
    tracers = required_biogeochemical_tracers(model)

    # if no parameters specified, get all of them (except `sinking_velocities`)
    if isnothing(parameters)
        parameters = filter(p -> p ∉ (:sinking_velocities,), fieldnames(biogeochemistry))
    end
    # retrieve parameter values from model struct
    p = [getfield(model, f) for f in fieldnames(biogeochemistry) if f ∈ parameters]

    function modelODE!(du, u, p, t)
        params = NamedTuple{parameters}(p)
        model = biogeochemistry(; params...)

        PAR = PAR_f(t)
        for (i, tracer) in enumerate(tracers)
            du[i] = model(Val(tracer), 0, 0, 0, t, u..., PAR)
        end
    end

    # make sure the initial values are passed in the right order
    u0 = [eval(:($init_conditions.$t)) for t in tracers]

    prob = ODEProblem(modelODE!, u0, tspan, p)

    return prob
end
