using OrdinaryDiffEq

using Oceananigans.Biogeochemistry:
    required_biogeochemical_tracers, required_biogeochemical_auxiliary_fields

"""
    bgc_to_ode(biogeochemistry, PAR_f, init_conditions, tspan, p=nothing)


"""
function bgc_to_ode(biogeochemistry, PAR_f, init_conditions, tspan, parameters=nothing)

    # retrieve tracer order
    model = biogeochemistry()
    tracers = required_biogeochemical_tracers(model)

    # if no parameters specified, get all of them
    if isnothing(parameters)
        parameters = fieldnames(biogeochemistry)
    end
    # retrive parameter values from model struct
    p = [getfield(model, f) for f in fieldnames(biogeochemistry) if f ∈ parameters]

    function modelODE!(du, u, p, t)
        params = NamedTuple{parameters}(p)
        model = biogeochemistry(; params...)

        # TODO: what if there are additional auxiliary fields here?
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
