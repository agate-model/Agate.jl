using SciMLBase: ODEProblem

"""Return an ODEProblem that evaluates an Agate BGC in a well-mixed box."""
function ode_problem(
    bgc,
    u0,
    tspan;
    p=nothing,
    active_parameters=nothing,
    auxiliary=(;),
    coordinates=(0, 0, 0),
)
    layout = active_parameter_layout(active_parameters)
    if !isempty(layout) && p === nothing
        throw(ArgumentError("`p` must be provided when `active_parameters` is non-empty."))
    end

    aux_names = required_biogeochemical_auxiliary_fields(bgc)
    tracer_names = required_biogeochemical_tracers(bgc)
    active_map = isempty(layout) ? nothing : active_parameter_map(bgc, layout)

    function rhs!(du, u, parameters_vector, t)
        parameters = active_map === nothing ? bgc.parameters : ActiveParameters(bgc.parameters, parameters_vector, active_map)
        aux_values = ode_auxiliary_values(aux_names, auxiliary, t)
        x, y, z = coordinates

        for (n, tracer) in enumerate(tracer_names)
            du[n] = evaluate_tendency(bgc, parameters, Val(tracer), x, y, z, t, u..., aux_values...)
        end

        return nothing
    end

    return p === nothing ? ODEProblem(rhs!, u0, tspan) : ODEProblem(rhs!, u0, tspan, p)
end

function ode_auxiliary_values(aux_names, auxiliary::NamedTuple, t)
    return map(aux_names) do name
        value = getproperty(auxiliary, name)
        return value isa Function ? value(t) : value
    end
end

ode_auxiliary_values(aux_names, auxiliary::Tuple, t) = auxiliary
