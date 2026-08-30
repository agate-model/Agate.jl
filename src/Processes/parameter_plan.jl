"""One model parameter after layout-dependent realization."""
struct PlannedParameter{D,A,S,L,M}
    name::Symbol
    definition::D
    rank::Int
    axes::A
    storage_shape::S
    storage_labels::L
    storage_diameters::M
    runtime_bound::Bool
end

"""Single host-side realization of parameter storage and runtime eligibility."""
struct ParameterPlan{P}
    parameters::P
end

_parameter_axis_entities(bindings::Tuple, binding_entities::Tuple, name::Symbol) = Tuple(
    entities for (binding, entities) in zip(bindings, binding_entities)
    if binding.parameter === name
)

function _layout_axis_entities(layout::ModelLayout, components::Tuple)
    entities = Symbol[]
    for component in components
        hasproperty(layout.component_entities, component) || throw(ArgumentError(
            "parameter applicability references unrealized component :$component",
        ))
        append!(entities, component_entities(layout, component))
    end
    return Tuple(entities)
end

function _resolved_binding_axis_entities(
    definition::CanonicalModelDefinition, layout::ModelLayout
)
    return map(definition.parameter_bindings) do binding
        map(components -> _layout_axis_entities(layout, components), binding.axis_components)
    end
end

function _layout_storage_order(layout::ModelLayout)
    labels = collect(layout.size_classes)
    seen = Set(labels)
    for entities in values(layout.component_entities), entity in entities
        entity in seen && continue
        push!(labels, entity)
        push!(seen, entity)
    end
    return Tuple(labels)
end

function _union_storage_labels(layout::ModelLayout, name::Symbol, rank::Int, axis_entities::Tuple)
    rank == 0 && return ()
    isempty(axis_entities) && throw(ArgumentError(
        "parameter :$name has slot-derived storage but no resolved process applicability",
    ))
    order = _layout_storage_order(layout)
    return ntuple(rank) do dimension
        flat = Set(entity for entities in axis_entities for entity in entities[dimension])
        labels = Tuple(label for label in order if label in flat)
        length(labels) == length(flat) || throw(ArgumentError(
            "parameter :$name applicability contains entities outside the realized layout",
        ))
        labels
    end
end

_planned_parameter_axes(definition, name, parameter::ConstructionParameter) =
    isnothing(parameter.axes) ? () : (parameter.axes,)

function _planned_parameter_axes(definition, name, parameter::Parameter)
    index = findfirst(binding -> binding.parameter === name, definition.parameter_bindings)
    isnothing(index) && throw(ArgumentError(
        "Parameter :$name has no scientific slot binding after canonicalization",
    ))
    return definition.parameter_bindings[index].axes
end

function _parameter_storage_labels(layout, name, axes, axis_entities)
    rank = length(axes)
    rank == 0 && return ()
    if isempty(axis_entities)
        axes == (:plankton,) || throw(ArgumentError(
            "construction-only parameter :$name has unsupported explicit axes $axes",
        ))
        return (layout.size_classes,)
    end
    return _union_storage_labels(layout, name, rank, axis_entities)
end

function _diameter_by_entity(layout::ModelLayout)
    values = Dict{Symbol,Any}(
        entity => (isfinite(diameter) && diameter > zero(diameter) ? diameter : nothing)
        for (entity, diameter) in zip(layout.size_classes, layout.size_class_diameters)
    )
    for component in keys(layout.component_entities)
        diameters = getproperty(layout.component_diameters, component)
        isnothing(diameters) && continue
        for (entity, diameter) in zip(getproperty(layout.component_entities, component), diameters)
            values[entity] = diameter
        end
    end
    return values
end

function _storage_diameters(rank, labels, diameter_by_entity)
    rank == 1 || return nothing
    axis = only(labels)
    all(label -> haskey(diameter_by_entity, label), axis) || return nothing
    return Tuple(diameter_by_entity[label] for label in axis)
end


function _planned_parameter(definition, layout, name, parameter, binding_entities, diameters)
    axis_entities = _parameter_axis_entities(
        definition.parameter_bindings, binding_entities, name
    )
    axes = _planned_parameter_axes(definition, name, parameter)
    rank = length(axes)
    labels = _parameter_storage_labels(layout, name, axes, axis_entities)
    return PlannedParameter(
        name,
        parameter,
        rank,
        axes,
        map(length, labels),
        labels,
        _storage_diameters(rank, labels, diameters),
        parameter isa Parameter,
    )
end

"""Build the authoritative realized parameter plan after `ModelLayout` exists."""
function build_parameter_plan(definition::CanonicalModelDefinition, layout::ModelLayout)
    definitions = isnothing(definition.parameters) ? (;) : definition.parameters
    binding_entities = _resolved_binding_axis_entities(definition, layout)
    diameters = _diameter_by_entity(layout)
    names = keys(definitions)
    parameters = NamedTuple{names}(ntuple(length(names)) do i
        name = names[i]
        _planned_parameter(
            definition, layout, name, getproperty(definitions, name), binding_entities, diameters
        )
    end)
    return ParameterPlan(parameters)
end

planned_parameter(plan::ParameterPlan, name::Symbol) =
    hasproperty(plan.parameters, name) ? getproperty(plan.parameters, name) :
    throw(ArgumentError("unknown parameter :$name"))

"""Resolve one realized entity directly to its realized parameter-storage index."""
function parameter_storage_index(
    plan::ParameterPlan, name::Symbol, dimension::Int, entity::Symbol
)
    parameter = planned_parameter(plan, name)
    1 <= dimension <= parameter.rank || throw(ArgumentError(
        "parameter :$name has no realized storage dimension $dimension",
    ))
    labels = parameter.storage_labels[dimension]
    index = findfirst(==(entity), labels)
    isnothing(index) && throw(ArgumentError(
        "parameter :$name realized entity :$entity is not present on realized storage dimension $dimension",
    ))
    return index
end

_runtime_parameter_names(plan::ParameterPlan) = Tuple(
    name for (name, parameter) in pairs(plan.parameters) if parameter.runtime_bound
)

"""Return the fully realized values that compiled runtime slots actually read."""
function runtime_parameter_values(plan::ParameterPlan, values::NamedTuple)
    names = _runtime_parameter_names(plan)
    return NamedTuple{names}(Tuple(getproperty(values, name) for name in names))
end

"""Compact host metadata used by introspection and active-parameter selection."""
function parameter_plan_metadata(plan::ParameterPlan)
    names = keys(plan.parameters)
    runtime_parameter_names = _runtime_parameter_names(plan)
    return NamedTuple{names}(ntuple(length(names)) do i
        parameter = getproperty(plan.parameters, names[i])
        derived_runtime_parameters = Tuple(
            target for target in runtime_parameter_names
            if begin
                default = getproperty(plan.parameters, target).definition.default
                default isa DerivedDefault && names[i] in default.deps
            end
        )
        (;
            rank=parameter.rank,
            axes=parameter.axes,
            shape=parameter.storage_shape,
            labels=parameter.storage_labels,
            runtime_bound=parameter.runtime_bound,
            derived_runtime_parameters,
        )
    end)
end
