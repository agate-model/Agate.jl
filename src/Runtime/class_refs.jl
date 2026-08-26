"""Reference one concrete class within a named population subgroup or component."""
struct ClassRef
    group::Symbol
    i::Int
end

"""Construct a `ClassRef` for host-side class lookup."""
@inline class(group::Symbol, i::Integer) = ClassRef(group, Int(i))

"""Return the number of classes in a realized component or population subgroup."""
function class_count(layout::ModelLayout, identity::Symbol)::Int
    if hasproperty(layout.component_classes, identity)
        return component_class_count(layout, identity)
    elseif hasproperty(layout.group_indices, identity)
        return length(getproperty(layout.group_indices, identity))
    end
    throw(ArgumentError("Unknown component or group symbol $identity"))
end

"""Resolve a `ClassRef` to one physical tracer position when that lookup is unambiguous."""
function resolve_class(layout::ModelLayout, cref::ClassRef)::Int
    if hasproperty(layout.component_classes, cref.group)
        nclasses = component_class_count(layout, cref.group)
        1 <= cref.i <= nclasses || throw(ArgumentError(
            "Class ordinal $(cref.i) is out of bounds for component $(cref.group) (valid 1:$nclasses).",
        ))
        state_mapping = component_state_tracers(layout, cref.group)
        if state_mapping === nothing
            return component_indices(layout, cref.group)[cref.i]
        end
        length(state_mapping) == 1 || throw(ArgumentError(
            "ClassRef for multi-state component $(cref.group) is ambiguous; select a prognostic state explicitly",
        ))
        return state_indices(layout, cref.group, only(keys(state_mapping)))[cref.i]
    elseif hasproperty(layout.group_indices, cref.group)
        global_indices = getproperty(layout.group_indices, cref.group)
        1 <= cref.i <= length(global_indices) || throw(ArgumentError(
            "Class ordinal $(cref.i) is out of bounds for group $(cref.group) (valid 1:$(length(global_indices))).",
        ))
        global_index = global_indices[cref.i]
        class = layout.class_symbols[global_index]
        hasproperty(layout.tracer_indices, class) || throw(ArgumentError(
            "ClassRef for group $(cref.group) is ambiguous because its population has multiple prognostic states",
        ))
        return getproperty(layout.tracer_indices, class)
    end
    throw(ArgumentError("Unknown component or group symbol $(cref.group)"))
end

"""Resolve one explicit population-state `ClassRef` to its physical tracer position."""
function resolve_state(layout::ModelLayout, cref::ClassRef, reference::PopulationStateRef)::Int
    cref.group === reference.population || throw(ArgumentError(
        "ClassRef component :$(cref.group) does not match state reference population :$(reference.population)",
    ))
    indices = state_indices(layout, reference.population, reference.state)
    1 <= cref.i <= length(indices) || throw(ArgumentError(
        "Class ordinal $(cref.i) is out of bounds for component $(cref.group) state $(reference.state) (valid 1:$(length(indices))).",
    ))
    return indices[cref.i]
end
