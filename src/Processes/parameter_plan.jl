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
        any(
            binding -> binding.parameter === name && binding.runtime_bound,
            definition.parameter_bindings,
        ),
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
    runtime_names = _runtime_parameter_names(plan)
    return NamedTuple{names}(ntuple(length(names)) do i
        parameter = getproperty(plan.parameters, names[i])
        derived_runtime_parameters = Tuple(
            target for target in runtime_names
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

function _realized_axis_parameter_value(
    plan::ParameterPlan,
    parameter_values::NamedTuple,
    binding::ParameterBinding,
    class::Symbol,
)
    length(binding.axes) == 1 || throw(ArgumentError(
        "parameter :$(binding.parameter) must have exactly one ecological storage axis for class-local validation",
    ))
    index = parameter_storage_index(plan, binding.parameter, 1, class)
    return getproperty(parameter_values, binding.parameter)[index]
end

function _validate_quota_bounds(
    definition::CanonicalModelDefinition,
    layout::ModelLayout,
    plan::ParameterPlan,
    parameter_values::NamedTuple,
    named::NamedProcess,
    path::Tuple,
    slot_refs::NamedTuple,
    population::Symbol,
)
    slots = _resolved_slot_bindings(definition, slot_refs)
    minimum_binding = slots.minimum_quota
    maximum_binding = slots.maximum_quota
    for class in component_classes(layout, population)
        minimum = _realized_axis_parameter_value(plan, parameter_values, minimum_binding, class)
        maximum = _realized_axis_parameter_value(plan, parameter_values, maximum_binding, class)
        context = "process :$(process_id(named)) path $path ecological class :$class"
        minimum > zero(minimum) || throw(ArgumentError(
            "$context parameter :$(minimum_binding.parameter) must be > 0; got $minimum",
        ))
        maximum > minimum || throw(ArgumentError(
            "$context parameter :$(maximum_binding.parameter) must be greater than :$(minimum_binding.parameter)=$minimum; got $maximum",
        ))
    end
    return nothing
end

function _validate_parameter_constraint(
    definition::CanonicalModelDefinition,
    layout::ModelLayout,
    plan::ParameterPlan,
    parameter_values::NamedTuple,
    named::NamedProcess,
    path::Tuple,
    slot_refs::NamedTuple,
    population::Symbol,
    slot::Symbol,
    rule::Symbol,
)
    binding = getproperty(_resolved_slot_bindings(definition, slot_refs), slot)
    for class in component_classes(layout, population)
        value = _realized_axis_parameter_value(plan, parameter_values, binding, class)
        valid = rule === :positive ? value > zero(value) : value >= zero(value)
        valid || throw(ArgumentError(
            "process :$(process_id(named)) path $path ecological class :$class " *
            "parameter :$(binding.parameter) must be $(rule === :positive ? "> 0" : ">= 0"); got $value",
        ))
    end
    return nothing
end

function _validate_quota_factor_science(
    definition, layout, plan, parameter_values, named, path, factor, refs
)
    if factor isa QuotaResponse
        _validate_quota_bounds(
            definition, layout, plan, parameter_values, named, path, refs.slots,
            factor.target.population,
        )
    end
    for (name, child) in pairs(factor_children(factor))
        _validate_quota_factor_science(
            definition,
            layout,
            plan,
            parameter_values,
            named,
            factor_child_path(path, factor, name),
            child,
            getproperty(refs.children, name),
        )
    end
    return nothing
end

function _validate_quota_science(
    definition::CanonicalModelDefinition,
    layout::ModelLayout,
    plan::ParameterPlan,
    parameter_values::NamedTuple,
)
    for named in values(definition.processes)
        for (name, factor) in pairs(factors(named))
            _validate_quota_factor_science(
                definition, layout, plan, parameter_values, named, (:factors, name), factor,
                getproperty(named.binding_refs.factors, name),
            )
        end
        process = named.process
        if process isa NutrientUptake
            path = ()
            _validate_quota_bounds(
                definition, layout, plan, parameter_values, named, path,
                named.binding_refs.process, process.population,
            )
            for (slot, rule) in (
                (:maximum_rate, :nonnegative), (:K, :nonnegative), (:hill, :positive),
            )
                _validate_parameter_constraint(
                    definition,
                    layout,
                    plan,
                    parameter_values,
                    named,
                    path,
                    named.binding_refs.process,
                    process.population,
                    slot,
                    rule,
                )
            end
        end
    end
    return nothing
end

function _validate_product_fractions(
    definition::CanonicalModelDefinition, parameter_values::NamedTuple
)
    for named in values(definition.processes)
        products = process_products(named.process)
        (isnothing(products) || length(products.targets) == 1) && continue
        names = keys(products.fractions)
        fractions = NamedTuple{names}(Tuple(begin
            refs = getproperty(named.binding_refs.products.fractions, product)
            binding = _resolved_slot_bindings(definition, refs).fraction
            value = getproperty(parameter_values, binding.parameter)
            value isa Real || throw(ArgumentError(
                "product fraction parameter :$(binding.parameter) for process :$(process_id(named)) must resolve to a scalar Real",
            ))
            zero(value) <= value <= one(value) || throw(ArgumentError(
                "product fraction parameter :$(binding.parameter) for process :$(process_id(named)) must lie in [0, 1]; got $value",
            ))
            value
        end for product in names))

        independent = Tuple(
            getproperty(fractions, product) for product in names if product !== products.balanced
        )
        total = sum(independent)
        balance = one(total) - total
        zero(balance) <= balance <= one(balance) || throw(ArgumentError(
            "product fractions for process :$(process_id(named)) leave conservative balance $balance for :$(products.balanced); expected a value in [0, 1]",
        ))

        if hasproperty(fractions, products.balanced)
            supplied = getproperty(fractions, products.balanced)
            tolerance = balance isa AbstractFloat ? 100 * eps(one(balance)) : zero(balance)
            isapprox(supplied, balance; rtol=zero(balance), atol=tolerance) || throw(ArgumentError(
                "explicit product fraction for :$(products.balanced) in process :$(process_id(named)) must agree with the conservative balance $balance; got $supplied",
            ))
        end
    end
    return nothing
end

"""Validate scientific constraints directly from canonical processes and realized values."""
function validate_realized_science(
    definition::CanonicalModelDefinition,
    layout::ModelLayout,
    plan::ParameterPlan,
    parameter_values::NamedTuple,
)
    _validate_product_fractions(definition, parameter_values)
    _validate_quota_science(definition, layout, plan, parameter_values)
    return nothing
end
