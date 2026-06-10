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
    active_map = active_parameters === nothing ? (;) : active_parameters.map
    if !isempty(active_map) && p === nothing
        throw(ArgumentError("`p` must be provided when `active_parameters` is non-empty."))
    end

    aux_names = required_biogeochemical_auxiliary_fields(bgc)
    aux_source = ode_auxiliary_source(aux_names, auxiliary)
    tracer_names = required_biogeochemical_tracers(bgc)

    if isempty(active_map)
        parameterized_bgc = ParameterizedBGC(bgc, bgc.parameters)

        rhs! = function (du, u, parameters_vector, t)
            aux_values = ode_auxiliary_values(aux_source, t)
            ode_tendencies!(du, parameterized_bgc, tracer_names, coordinates, u, aux_values, t)
            return nothing
        end
    else
        rhs! = function (du, u, parameters_vector, t)
            parameters = ActiveParameters(bgc.parameters, parameters_vector, active_map)
            parameterized_bgc = ParameterizedBGC(bgc, parameters)
            aux_values = ode_auxiliary_values(aux_source, t)
            ode_tendencies!(du, parameterized_bgc, tracer_names, coordinates, u, aux_values, t)
            return nothing
        end
    end

    return p === nothing ? ODEProblem(rhs!, u0, tspan) : ODEProblem(rhs!, u0, tspan, p)
end

@inline function ode_tendencies!(du, bgc, tracer_names, coordinates, u, aux_values, t)
    x, y, z = coordinates

    for (n, tracer) in enumerate(tracer_names)
        du[n] = evaluate_tendency(bgc, Val(tracer), x, y, z, t, u..., aux_values...)
    end

    return nothing
end

struct ODEAuxiliarySource{V,E}
    values::V
end

function ode_auxiliary_source(aux_names, auxiliary::NamedTuple)
    values = map(aux_names) do name
        hasproperty(auxiliary, name) || throw(ArgumentError("Missing auxiliary value :$name."))
        return getproperty(auxiliary, name)
    end
    return ODEAuxiliarySource{typeof(values), true}(values)
end

ode_auxiliary_source(aux_names, auxiliary::Tuple) = ODEAuxiliarySource{typeof(auxiliary), false}(auxiliary)

@inline function ode_auxiliary_values(source::ODEAuxiliarySource{V, true}, t) where {V}
    return map(source.values) do value
        value isa Function ? value(t) : value
    end
end

@inline ode_auxiliary_values(source::ODEAuxiliarySource{V, false}, t) where {V} = source.values
