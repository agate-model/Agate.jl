# Canonical process identity and drivers

function _canonicalize_process(id::Symbol, process, components::NamedTuple)
    process isa AbstractProcess || throw(
        ArgumentError("process :$id must be an AbstractProcess; got $(typeof(process))"),
    )
    component_names = keys(components)
    missing_names = unique(Tuple(
        reference for reference in _collect_process_component_references(process)
        if !(reference in component_names)
    ))
    isempty(missing_names) || throw(
        ArgumentError("process :$id references unknown components $missing_names"),
    )
    semantic_facts = process_facts(process, id, components)
    for (name, factor) in pairs(factors(process))
        _validate_factor_for_process(id, process, factor, components, (:factors, name))
    end
    return CanonicalProcess(id, process, semantic_facts)
end

function _canonical_processes(processes::NamedTuple, components::NamedTuple)
    names = sort!(collect(keys(processes)); by=String)
    names_tuple = Tuple(names)
    process_values = Tuple(
        _canonicalize_process(name, getproperty(processes, name), components) for name in names
    )
    return NamedTuple{names_tuple}(process_values)
end

function _collect_driver_identities!(identities::Vector{Symbol}, factor::AbstractFactor)
    for input in factor_inputs(factor)
        input isa FactorDriver || continue
        input.identity in identities || push!(identities, input.identity)
    end
    for subfactor in values(_resolve_factor_subfactors(factor))
        _collect_driver_identities!(identities, subfactor)
    end
    return nothing
end

function _canonical_driver_identities(processes::NamedTuple)
    identities = Symbol[]
    for process in values(processes), factor in values(factors(process))
        _collect_driver_identities!(identities, factor)
    end
    sort!(identities; by=String)
    return Tuple(identities)
end


# Parameter definition and binding resolution

function _resolve_binding_value(bindings::NamedTuple, slot::ParameterSlot, qualifier)
    explicit = hasproperty(bindings, slot.name)
    value = explicit ? getproperty(bindings, slot.name) : slot.name
    !isnothing(slot.qualify) && isnothing(qualifier) && throw(ArgumentError(
        "parameter slot :$(slot.name) requires qualifier :$(slot.qualify), but no qualifier context was resolved",
    ))
    value isa Symbol && return value, explicit, nothing

    value isa NamedTuple || throw(ArgumentError(
        "binding :$(slot.name) must be a parameter Symbol or one-level qualifier NamedTuple",
    ))
    isnothing(slot.qualify) && throw(ArgumentError(
        "binding :$(slot.name) uses a qualifier map but the slot is unqualified",
    ))
    qualifier isa Qualifier && qualifier.axis === slot.qualify || throw(ArgumentError(
        "binding :$(slot.name) requires qualifier :$(slot.qualify)",
    ))
    qualifier_value = qualifier.value
    hasproperty(value, qualifier_value) || throw(ArgumentError(
        "binding :$(slot.name) has no entry for qualifier :$qualifier_value",
    ))
    parameter = getproperty(value, qualifier_value)
    parameter isa Symbol || throw(ArgumentError(
        "binding :$(slot.name) qualifier :$qualifier_value must map to a parameter Symbol",
    ))
    return parameter, true, keys(value)
end

function _collect_parameter_uses(processes::NamedTuple)
    uses = Any[]
    seen = Set{Any}()
    names = keys(processes)
    refs = NamedTuple{names}(Tuple(
        _visit_process_slots!(uses, seen, named) for named in values(processes)
    ))
    return Tuple(uses), refs
end

function _resolve_parameter_definitions(definitions)
    isnothing(definitions) && return nothing, Set{Symbol}()
    definitions isa NamedTuple || throw(
        ArgumentError("model parameters must be a NamedTuple of Parameter or ConstructionParameter values"),
    )
    all(parameter -> parameter isa Union{Parameter,ConstructionParameter}, values(definitions)) || throw(
        ArgumentError("model parameters must contain only Parameter or ConstructionParameter values"),
    )

    names = Tuple(sort!(collect(keys(definitions)); by=String))
    definitions = NamedTuple{names}(
        Tuple(getproperty(definitions, name) for name in names)
    )
    definition_set = Set(names)

    dependency_names = Set{Symbol}()
    for (name, parameter) in pairs(definitions)
        default = parameter.default
        default isa DerivedDefault || continue
        for dependency in default.deps
            dependency in definition_set || throw(
                ArgumentError(
                    "parameter :$name default depends on undeclared parameter :$dependency",
                ),
            )
            dependency_default = getproperty(definitions, dependency).default
            dependency_default isa DerivedDefault && throw(
                ArgumentError(
                    "parameter :$name default depends on derived parameter :$dependency; " *
                    "DerivedDefault dependencies must be non-derived parameters",
                ),
            )
            push!(dependency_names, dependency)
        end
    end
    return definitions, dependency_names
end

_parameter_binding(use) = ParameterBinding(
    use.process, use.path, use.slot, use.axes, use.axis_components, use.parameter, use.domain
)

function _resolve_parameter_bindings(
    uses::Tuple, definitions, dependency_names::Set{Symbol}
)
    for use in uses
        isnothing(use.binding_qualifiers) && continue
        consumed = Tuple(other.qualifier.value for other in uses if
            other.process === use.process && other.path == use.path &&
            other.slot === use.slot && other.qualifier isa Qualifier)
        extra = Tuple(key for key in use.binding_qualifiers if !(key in consumed))
        isempty(extra) || throw(ArgumentError("binding :$(use.slot) has unused qualifier entries $extra"))
    end

    if isnothing(definitions)
        bindings = Tuple(_parameter_binding(use) for use in uses)
        return bindings
    end

    definition_names = Set(keys(definitions))
    for use in uses
        use.parameter in definition_names || throw(ArgumentError(
            "inline binding for process :$(use.process) path $(use.path) slot :$(use.slot) " *
            "qualifier $(use.qualifier) names undeclared parameter :$(use.parameter)",
        ))
    end

    by_parameter = Dict{Symbol,Vector{Any}}()
    for use in uses
        push!(get!(by_parameter, use.parameter, Any[]), use)
    end
    for (parameter, parameter_uses) in by_parameter
        if length(parameter_uses) > 1 && any(use -> !use.explicit, parameter_uses)
            locations = Tuple(
                (process=use.process, path=use.path, slot=use.slot, qualifier=use.qualifier)
                for use in parameter_uses
            )
            throw(ArgumentError(
                "parameter :$parameter is implicitly bound by multiple parameter slots $locations; " *
                "bind every shared use explicitly",
            ))
        end
        required_axes = unique(Tuple(use.axes for use in parameter_uses))
        length(required_axes) == 1 || throw(ArgumentError(
            "parameter :$parameter is bound to slots with incompatible semantic axes $(Tuple(required_axes))",
        ))
        required_domains = unique(Tuple(use.domain for use in parameter_uses))
        length(required_domains) == 1 || throw(ArgumentError(
            "parameter :$parameter is bound to slots with incompatible scientific domains $(Tuple(required_domains))",
        ))
    end

    for (name, parameter) in pairs(definitions)
        parameter_uses = get(by_parameter, name, Any[])
        if parameter isa ConstructionParameter
            isempty(parameter_uses) || throw(ArgumentError(
                "ConstructionParameter :$name is construction-only and cannot bind to a process slot; use Parameter for runtime process parameters",
            ))
            name in dependency_names || throw(ArgumentError(
                "ConstructionParameter :$name must be used by at least one DerivedDefault",
            ))
        else
            isempty(parameter_uses) && throw(ArgumentError(
                "Parameter :$name is not bound to a process slot; use ConstructionParameter for construction-only DerivedDefault inputs",
            ))
        end
    end

    return Tuple(_parameter_binding(use) for use in uses)
end

function _attach_binding_refs(processes::NamedTuple, refs::NamedTuple)
    names = keys(processes)
    return NamedTuple{names}(Tuple(
        begin
            named = getproperty(processes, name)
            CanonicalProcess(named.id, named.process, named.semantic_facts, getproperty(refs, name))
        end
        for name in names
    ))
end

# Complete model canonicalization

"""Canonicalize process identity and resolve inline parameter bindings.

Process instances and parameter definitions are canonicalized by stable identity, so their
declaration order does not change canonical scientific identity. Component order is preserved
within component kinds because `ModelLayout` exposes concrete runtime tracer order. Local
formulation slots bind directly to stable model parameter names during canonicalization.
"""
function canonicalize_model(definition::ModelDefinition)
    all(component -> component isa Union{Plankton,Pool}, values(definition.components)) ||
        throw(ArgumentError("model components must be Plankton or Pool values"))
    canonical_processes = _canonical_processes(
        definition.processes, definition.components
    )
    parameters, dependency_names = _resolve_parameter_definitions(definition.parameters)
    uses, binding_refs = _collect_parameter_uses(canonical_processes)
    bindings = _resolve_parameter_bindings(uses, parameters, dependency_names)
    bound_processes = _attach_binding_refs(canonical_processes, binding_refs)
    return CanonicalModelDefinition(
        definition.components,
        bound_processes,
        parameters,
        _canonical_driver_identities(canonical_processes),
        bindings,
    )
end
