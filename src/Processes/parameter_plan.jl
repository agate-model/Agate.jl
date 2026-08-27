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

"""Realized product-fraction constraint evaluated after parameter materialization."""
struct ProductFractionCheck{F}
    process::Symbol
    balanced::Symbol
    fractions::F
end

"""One class-local minimum/maximum quota inequality check."""
struct QuotaBoundsCheck{I,J}
    process::Symbol
    path::Tuple
    class::Symbol
    minimum_parameter::Symbol
    minimum_index::I
    maximum_parameter::Symbol
    maximum_index::J
end

"""One class-local scalar inequality for a realized model parameter."""
struct ParameterConstraintCheck{I}
    process::Symbol
    path::Tuple
    class::Symbol
    parameter::Symbol
    index::I
    rule::Symbol
end

"""Single host-side realization of parameter storage, indexing, and runtime eligibility."""
struct ParameterPlan{P,L,R,C}
    parameters::P
    slot_lookup::L
    runtime_names::R
    science_checks::C
end

_parameter_axis_classes(bindings::Tuple, binding_classes::Tuple, name::Symbol) = Tuple(
    classes for (binding, classes) in zip(bindings, binding_classes)
    if binding.parameter === name
)

function _binding_axis_components(
    process::NamedProcess, binding::ParameterBinding, axis::Symbol
)
    qualifier = binding.qualifier
    qualifier isa Qualifier && qualifier.axis === axis && return (qualifier.value,)
    process_participants = participants(process)
    hasproperty(process_participants, axis) || throw(ArgumentError(
        "parameter applicability axis :$axis is not a participant role of process :$(process_id(process))",
    ))
    return getproperty(process_participants, axis)
end

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

function _binding_applicability_axes(process::NamedProcess, binding::ParameterBinding)
    isempty(binding.axes) || return binding.axes
    qualifier = binding.qualifier
    qualifier isa Qualifier || return ()
    return hasproperty(participants(process), qualifier.axis) ? (qualifier.axis,) : ()
end

function _resolved_binding_axis_classes(
    definition::NormalizedModelDefinition, layout::ModelLayout
)
    return map(definition.parameter_bindings) do binding
        process = getproperty(definition.processes, binding.process)
        axes = _binding_applicability_axes(process, binding)
        map(axis -> _layout_axis_classes(layout, _binding_axis_components(process, binding, axis)), axes)
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

_explicit_axis_labels(layout::ModelLayout, axis::Symbol) =
    axis === :plankton ? layout.class_symbols :
    Tuple(layout.class_symbols[index] for index in axis_indices(layout, axis))

_axis_tuple(::Nothing) = ()
_axis_tuple(axis::Symbol) = (axis,)
_axis_tuple(axes::Tuple) = axes

_planned_parameter_axes(definition, name, parameter::MetaParameter) = _axis_tuple(parameter.axes)

function _planned_parameter_axes(definition, name, parameter::Parameter)
    index = findfirst(binding -> binding.parameter === name, definition.parameter_bindings)
    isnothing(index) && throw(ArgumentError(
        "Parameter :$name has no scientific slot binding after normalization",
    ))
    return definition.parameter_bindings[index].axes
end

function _parameter_storage_labels(layout, name, axes, axis_classes)
    rank = length(axes)
    rank == 0 && return ()
    isempty(axis_classes) && return map(axis -> _explicit_axis_labels(layout, axis), axes)
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


function _runtime_binding(definition::NormalizedModelDefinition, binding::ParameterBinding)
    binding.slot === :fraction || return true
    named = getproperty(definition.processes, binding.process)
    products = process_products(named.process)
    products isa Products || return true
    binding.path == product_path(named.process) || return true
    qualifier = binding.qualifier
    return !(qualifier isa Qualifier && qualifier.value === products.balanced)
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
            binding -> binding.parameter === name && _runtime_binding(definition, binding),
            definition.parameter_bindings,
        ),
    )
end

function _storage_index_map(name, local_labels, storage_labels)
    return Tuple(map(local_labels) do label
        index = findfirst(==(label), storage_labels)
        isnothing(index) && throw(ArgumentError(
            "parameter :$name process class :$label is not present on its realized storage axis",
        ))
        index
    end)
end

function _slot_index_maps(binding, axis_classes, parameter)
    return ntuple(length(binding.axes)) do dimension
        _storage_index_map(
            binding.parameter,
            axis_classes[dimension],
            parameter.storage_labels[dimension],
        )
    end
end


function _slot_axis_index_map(slot_lookup, binding::ParameterBinding)
    mappings = get(slot_lookup, _binding_key(binding)) do
        throw(ArgumentError(
            "parameter plan has no realized slot for :$(binding.parameter) at process :$(binding.process)",
        ))
    end
    length(mappings) == 1 || throw(ArgumentError(
        "quota parameter :$(binding.parameter) must have exactly one population axis",
    ))
    return only(mappings)
end

function _append_quota_bounds_checks!(
    checks, definition, layout, slot_lookup, named, path, node, population::Symbol
)
    slots = parameter_slot_bindings(definition, named, path, node)
    minimum_map = _slot_axis_index_map(slot_lookup, slots.minimum_quota)
    maximum_map = _slot_axis_index_map(slot_lookup, slots.maximum_quota)
    classes = component_classes(layout, population)
    length(classes) == length(minimum_map) == length(maximum_map) || throw(ArgumentError(
        "quota parameter realization for process :$(process_id(named)) does not match population :$population classes",
    ))
    for (local_index, class) in pairs(classes)
        push!(checks, QuotaBoundsCheck(
            process_id(named), path, class,
            slots.minimum_quota.parameter, (minimum_map[local_index],),
            slots.maximum_quota.parameter, (maximum_map[local_index],),
        ))
    end
    return nothing
end

function _append_parameter_constraint_checks!(
    checks, definition, layout, slot_lookup, named, path, node,
    population::Symbol, slot::Symbol, rule::Symbol,
)
    binding = getproperty(parameter_slot_bindings(definition, named, path, node), slot)
    mapping = _slot_axis_index_map(slot_lookup, binding)
    classes = component_classes(layout, population)
    length(classes) == length(mapping) || throw(ArgumentError(
        "parameter :$(binding.parameter) realization for process :$(process_id(named)) does not match population :$population classes",
    ))
    for (local_index, class) in pairs(classes)
        push!(checks, ParameterConstraintCheck(
            process_id(named), path, class, binding.parameter, (mapping[local_index],), rule,
        ))
    end
    return nothing
end

function _quota_factor_checks!(checks, definition, layout, slot_lookup, named, path, factor)
    if factor isa QuotaResponse
        _append_quota_bounds_checks!(
            checks, definition, layout, slot_lookup, named, path, factor, factor.target.population
        )
    end
    for (name, child) in pairs(factor_children(factor))
        _quota_factor_checks!(
            checks, definition, layout, slot_lookup, named,
            factor_child_path(path, factor, name), child,
        )
    end
    return nothing
end

function _quota_science_checks(definition, layout, slot_lookup)
    checks = Any[]
    for named in values(definition.processes)
        for (name, factor) in pairs(factors(named))
            _quota_factor_checks!(
                checks, definition, layout, slot_lookup, named, (:factors, name), factor
            )
        end
        process = named.process
        if process isa NutrientUptake
            path = ()
            _append_quota_bounds_checks!(
                checks, definition, layout, slot_lookup, named, path, process, process.population
            )
            for (slot, rule) in (
                (:maximum_rate, :nonnegative), (:K, :nonnegative), (:hill, :positive),
            )
                _append_parameter_constraint_checks!(
                    checks, definition, layout, slot_lookup, named, path, process,
                    process.population, slot, rule,
                )
            end
        end
    end
    return Tuple(checks)
end

function _product_fraction_checks(definition::NormalizedModelDefinition)
    checks = Any[]
    for named in values(definition.processes)
        products = process_products(named.process)
        (isnothing(products) || length(products.targets) == 1) && continue
        names = keys(products.fractions)
        fractions = NamedTuple{names}(Tuple(
            parameter_slot_bindings(
                definition, named, product_path(named.process), products;
                context=(product=product,),
            ).fraction.parameter
            for product in names
        ))
        push!(checks, ProductFractionCheck(process_id(named), products.balanced, fractions))
    end
    return Tuple(checks)
end

"""Build the authoritative realized parameter plan after `ModelLayout` exists."""
function build_parameter_plan(definition::NormalizedModelDefinition, layout::ModelLayout)
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
    slot_lookup = Dict(
        _binding_key(binding) => _slot_index_maps(
            binding, classes, getproperty(parameters, binding.parameter)
        )
        for (binding, classes) in zip(definition.parameter_bindings, binding_classes)
    )
    runtime_names = Tuple(name for name in names if getproperty(parameters, name).runtime_bound)
    science_checks = (
        _product_fraction_checks(definition)...,
        _quota_science_checks(definition, layout, slot_lookup)...,
    )
    return ParameterPlan(parameters, slot_lookup, runtime_names, science_checks)
end

planned_parameter(plan::ParameterPlan, name::Symbol) =
    hasproperty(plan.parameters, name) ? getproperty(plan.parameters, name) :
    throw(ArgumentError("unknown parameter :$name"))

function planned_parameter_slot(plan::ParameterPlan, binding::ParameterBinding)
    return get(plan.slot_lookup, _binding_key(binding)) do
        throw(ArgumentError(
            "parameter plan has no realized slot for :$(binding.parameter) at process :$(binding.process)",
        ))
    end
end

"""Return the fully realized values that compiled runtime slots actually read."""
function runtime_parameter_values(plan::ParameterPlan, values::NamedTuple)
    return NamedTuple{plan.runtime_names}(
        Tuple(getproperty(values, name) for name in plan.runtime_names)
    )
end

"""Compact host metadata used by introspection and active-parameter selection."""
function parameter_plan_metadata(plan::ParameterPlan)
    names = keys(plan.parameters)
    return NamedTuple{names}(ntuple(length(names)) do i
        parameter = getproperty(plan.parameters, names[i])
        derived_runtime_parameters = Tuple(
            target for target in plan.runtime_names
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


@inline function _realized_parameter_value(values::NamedTuple, parameter::Symbol, index::Tuple)
    value = getproperty(values, parameter)
    return isempty(index) ? value : value[index...]
end

function validate_realized_science(check::QuotaBoundsCheck, values::NamedTuple)
    minimum = _realized_parameter_value(values, check.minimum_parameter, check.minimum_index)
    maximum = _realized_parameter_value(values, check.maximum_parameter, check.maximum_index)
    context = "process :$(check.process) path $(check.path) ecological class :$(check.class)"
    minimum > zero(minimum) || throw(ArgumentError(
        "$context parameter :$(check.minimum_parameter) must be > 0; got $minimum",
    ))
    maximum > minimum || throw(ArgumentError(
        "$context parameter :$(check.maximum_parameter) must be greater than :$(check.minimum_parameter)=$minimum; got $maximum",
    ))
    return nothing
end

function validate_realized_science(check::ParameterConstraintCheck, values::NamedTuple)
    value = _realized_parameter_value(values, check.parameter, check.index)
    valid = check.rule === :positive ? value > zero(value) : value >= zero(value)
    valid || throw(ArgumentError(
        "process :$(check.process) path $(check.path) ecological class :$(check.class) " *
        "parameter :$(check.parameter) must be $(check.rule === :positive ? "> 0" : ">= 0"); got $value",
    ))
    return nothing
end

function validate_realized_science(check::ProductFractionCheck, values::NamedTuple)
    names = keys(check.fractions)
    fractions = NamedTuple{names}(Tuple(begin
        parameter = getproperty(check.fractions, product)
        value = getproperty(values, parameter)
        value isa Real || throw(ArgumentError(
            "product fraction parameter :$parameter for process :$(check.process) must resolve to a scalar Real",
        ))
        zero(value) <= value <= one(value) || throw(ArgumentError(
            "product fraction parameter :$parameter for process :$(check.process) must lie in [0, 1]; got $value",
        ))
        value
    end for product in names))

    independent = Tuple(getproperty(fractions, product) for product in names if product !== check.balanced)
    total = sum(independent)
    balance = one(total) - total
    zero(balance) <= balance <= one(balance) || throw(ArgumentError(
        "product fractions for process :$(check.process) leave conservative balance $balance for :$(check.balanced); expected a value in [0, 1]",
    ))

    if hasproperty(fractions, check.balanced)
        supplied = getproperty(fractions, check.balanced)
        tolerance = balance isa AbstractFloat ? 100 * eps(one(balance)) : zero(balance)
        isapprox(supplied, balance; rtol=zero(balance), atol=tolerance) || throw(ArgumentError(
            "explicit product fraction for :$(check.balanced) in process :$(check.process) must agree with the conservative balance $balance; got $supplied",
        ))
    end
    return nothing
end

"""Validate realized scientific parameter constraints through one dispatched pass."""
function validate_realized_science(plan::ParameterPlan, values::NamedTuple)
    foreach(check -> validate_realized_science(check, values), plan.science_checks)
    return nothing
end
