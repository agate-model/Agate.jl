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

_parameter_axis_classes(bindings::Tuple, binding_classes::Tuple, name::Symbol) = Tuple(
    classes for (binding, classes) in zip(bindings, binding_classes)
    if binding.parameter === name
)

function _layout_axis_classes(layout::ModelLayout, components::Tuple)
    classes = Symbol[]
    for component in components
        hasproperty(layout.component_classes, component) || throw(ArgumentError(
            "parameter applicability references unrealized component :$component",
        ))
        append!(classes, component_classes(layout, component))
    end
    return Tuple(classes)
end

function _resolved_binding_axis_classes(
    definition::CanonicalModelDefinition, layout::ModelLayout
)
    return map(definition.parameter_bindings) do binding
        map(components -> _layout_axis_classes(layout, components), binding.axis_components)
    end
end

function _layout_storage_order(layout::ModelLayout)
    labels = collect(layout.class_symbols)
    seen = Set(labels)
    for classes in values(layout.component_classes), class in classes
        class in seen && continue
        push!(labels, class)
        push!(seen, class)
    end
    return Tuple(labels)
end

function _union_storage_labels(layout::ModelLayout, name::Symbol, rank::Int, axis_classes::Tuple)
    rank == 0 && return ()
    isempty(axis_classes) && throw(ArgumentError(
        "parameter :$name has slot-derived storage but no resolved process applicability",
    ))
    order = _layout_storage_order(layout)
    return ntuple(rank) do dimension
        flat = Set(class for classes in axis_classes for class in classes[dimension])
        labels = Tuple(label for label in order if label in flat)
        length(labels) == length(flat) || throw(ArgumentError(
            "parameter :$name applicability contains classes outside the realized layout",
        ))
        labels
    end
end

_planned_parameter_axes(definition, name, parameter::MetaParameter) =
    isnothing(parameter.axes) ? () : (parameter.axes,)

function _planned_parameter_axes(definition, name, parameter::Parameter)
    index = findfirst(binding -> binding.parameter === name, definition.parameter_bindings)
    isnothing(index) && throw(ArgumentError(
        "Parameter :$name has no scientific slot binding after canonicalization",
    ))
    return definition.parameter_bindings[index].axes
end

function _parameter_storage_labels(layout, name, axes, axis_classes)
    rank = length(axes)
    rank == 0 && return ()
    if isempty(axis_classes)
        axes == (:plankton,) || throw(ArgumentError(
            "construction-only parameter :$name has unsupported explicit axes $axes",
        ))
        return (layout.class_symbols,)
    end
    return _union_storage_labels(layout, name, rank, axis_classes)
end

function _diameter_by_class(layout::ModelLayout)
    values = Dict{Symbol,Any}(zip(layout.class_symbols, layout.diameters))
    for component in keys(layout.component_classes)
        diameters = getproperty(layout.component_diameters, component)
        isnothing(diameters) && continue
        for (class, diameter) in zip(getproperty(layout.component_classes, component), diameters)
            values[class] = diameter
        end
    end
    return values
end

function _storage_diameters(rank, labels, diameter_by_class)
    rank == 1 || return nothing
    axis = only(labels)
    all(label -> haskey(diameter_by_class, label), axis) || return nothing
    return Tuple(diameter_by_class[label] for label in axis)
end


function _planned_parameter(definition, layout, name, parameter, binding_classes, diameters)
    axis_classes = _parameter_axis_classes(
        definition.parameter_bindings, binding_classes, name
    )
    axes = _planned_parameter_axes(definition, name, parameter)
    rank = length(axes)
    labels = _parameter_storage_labels(layout, name, axes, axis_classes)
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
    binding_classes = _resolved_binding_axis_classes(definition, layout)
    diameters = _diameter_by_class(layout)
    names = keys(definitions)
    parameters = NamedTuple{names}(ntuple(length(names)) do i
        name = names[i]
        _planned_parameter(
            definition, layout, name, getproperty(definitions, name), binding_classes, diameters
        )
    end)
    return ParameterPlan(parameters)
end

planned_parameter(plan::ParameterPlan, name::Symbol) =
    hasproperty(plan.parameters, name) ? getproperty(plan.parameters, name) :
    throw(ArgumentError("unknown parameter :$name"))

"""Resolve one ecological class directly to its realized parameter-storage index."""
function parameter_storage_index(
    plan::ParameterPlan, name::Symbol, dimension::Int, class::Symbol
)
    parameter = planned_parameter(plan, name)
    1 <= dimension <= parameter.rank || throw(ArgumentError(
        "parameter :$name has no realized storage dimension $dimension",
    ))
    labels = parameter.storage_labels[dimension]
    index = findfirst(==(class), labels)
    isnothing(index) && throw(ArgumentError(
        "parameter :$name ecological class :$class is not present on realized storage dimension $dimension",
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
            shape=parameter.storage_shape,
            labels=parameter.storage_labels,
            runtime_bound=parameter.runtime_bound,
            derived_runtime_parameters,
        )
    end)
end

"""Return host metadata for runtime parameters stored on canonical consumer-by-prey axes."""
function interaction_axis_metadata(plan::ParameterPlan, layout::ModelLayout)
    consumers = Tuple(layout.class_symbols[index] for index in layout.consumer_indices)
    prey = Tuple(layout.class_symbols[index] for index in layout.prey_indices)
    names = Tuple(
        name for name in keys(plan.parameters)
        if begin
            parameter = getproperty(plan.parameters, name)
            parameter.runtime_bound &&
                parameter.axes == (:consumer, :resource) &&
                parameter.storage_labels == (consumers, prey)
        end
    )
    isempty(names) && return nothing
    return (; parameters=names, consumers, prey)
end
