using SciMLBase: ODEProblem

"""Return an ODEProblem that evaluates an Agate BGC in a well-mixed box."""
function ode_problem(
    bgc,
    u0,
    tspan;
    p=nothing,
    active_parameters=(;),
    auxiliary=(;),
    coordinates=(0, 0, 0),
)
    aux_names = required_biogeochemical_auxiliary_fields(bgc)

    function rhs!(du, u, parameters_vector, t)
        bgc_t = if isempty(active_parameters)
            bgc
        else
            parameterized(bgc, parameters_vector; active_parameters)
        end

        aux_values = ode_auxiliary_values(aux_names, auxiliary, t)
        x, y, z = coordinates
        tracer_names = required_biogeochemical_tracers(bgc)

        for (n, tracer) in enumerate(tracer_names)
            du[n] = bgc_t(Val(tracer), x, y, z, t, u..., aux_values...)
        end

        return nothing
    end

    if !isempty(active_parameters) && p === nothing
        throw(ArgumentError("`p` must be provided when `active_parameters` is non-empty."))
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
