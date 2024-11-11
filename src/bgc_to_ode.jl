using OrdinaryDiffEq

using Oceananigans.Biogeochemistry: required_biogeochemical_tracers

function bgc_to_ode(biogeochemistry, PAR_f, init_conditions, tspan, p=nothing)
    model = biogeochemistry()
    tracers = required_biogeochemical_tracers(model)

    function modelODE!(du, u, p, t)
        model = biogeochemistry(p...)

        PAR = PAR_f(t)
        for (i, tracer) in enumerate(tracers)
            du[i] = model(Val(tracer), 0, 0, 0, t, u..., PAR)
        end
    end

    u0 = [eval(:($init_conditions.$t)) for t in tracers]

    if isnothing(p)
        p = [getfield(model, f) for f in fieldnames(typeof(model))]
    end

    prob = ODEProblem(modelODE!, u0, tspan, p)

    return prob
end
