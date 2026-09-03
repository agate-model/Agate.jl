@inline function parameter_domain_valid(value, domain::Symbol)
    value isa Real && !(value isa Bool) && isfinite(value) || return false
    domain === :finite && return true
    domain === :nonnegative && return value >= zero(value)
    domain === :positive && return value > zero(value)
    domain === :unit_interval && return zero(value) <= value <= one(value)
    return false
end

_parameter_entries(value::Number) = (value,)
_parameter_entries(value) = value

function _validate_parameter_domains(
    definition::CanonicalModelDefinition, parameter_values::NamedTuple
)
    for binding in definition.parameter_bindings
        value = getproperty(parameter_values, binding.parameter)
        for entry in _parameter_entries(value)
            parameter_domain_valid(entry, binding.domain) && continue
            throw(ArgumentError(
                "process :$(binding.process) path $(binding.path) slot :$(binding.slot) " *
                "parameter :$(binding.parameter) must satisfy domain :$(binding.domain); got $entry",
            ))
        end
    end
    return nothing
end

"""One fully realized runtime parameter entry used by a coupled scientific constraint."""
struct ParameterValueRef{Parameter,Indices}
    indices::Indices
end

ParameterValueRef(parameter::Symbol, indices) = ParameterValueRef{parameter,typeof(indices)}(indices)
_parameter_name(::ParameterValueRef{Parameter}) where {Parameter} = Parameter

struct OrderedParameterConstraint{Lower,Upper,Path}
    lower::Lower
    upper::Upper
    process::Symbol
    path::Path
    entity::Symbol
end

struct ProductFractionConstraint{Terms}
    terms::Terms
    process::Symbol
end

@inline function _parameter_value(parameters, reference::ParameterValueRef{Parameter}) where {Parameter}
    value = getproperty(parameters, Parameter)
    isempty(reference.indices) && return value
    return value[reference.indices...]
end

function _validate_parameter_constraint(parameters, constraint::OrderedParameterConstraint)
    minimum = _parameter_value(parameters, constraint.lower)
    maximum = _parameter_value(parameters, constraint.upper)
    maximum > minimum && return nothing
    throw(ArgumentError(
        "process :$(constraint.process) path $(constraint.path) entity :$(constraint.entity) " *
        "parameter :$(_parameter_name(constraint.upper)) must be greater than " *
        ":$(_parameter_name(constraint.lower))=$minimum; got $maximum",
    ))
end

function _validate_parameter_constraint(parameters, constraint::ProductFractionConstraint)
    fractions = map(reference -> _parameter_value(parameters, reference), constraint.terms)
    total = sum(fractions)
    zero(total) <= total <= one(total) && return nothing
    throw(ArgumentError(
        "product fractions for process :$(constraint.process) sum to $total; expected a value " *
        "in [0, 1] so the omitted product receives a conservative remainder",
    ))
end

"""Validate coupled formulation constraints against any runtime-compatible parameter container."""
function validate_parameter_constraints(parameters, constraints::Tuple)
    for constraint in constraints
        _validate_parameter_constraint(parameters, constraint)
    end
    return nothing
end

function _append_quota_constraints!(
    constraints,
    definition::CanonicalModelDefinition,
    layout::ModelLayout,
    plan::ParameterPlan,
    named::CanonicalProcess,
    path::Tuple,
    slot_refs::NamedTuple,
)
    slots = _resolved_slot_bindings(definition, slot_refs)
    minimum_binding = slots.minimum_quota
    maximum_binding = slots.maximum_quota
    minimum_entities = _layout_axis_entities(layout, only(minimum_binding.axis_components))
    maximum_entities = _layout_axis_entities(layout, only(maximum_binding.axis_components))
    minimum_entities == maximum_entities || throw(ArgumentError(
        "quota bounds for process :$(process_id(named)) path $path resolve to different plankton entities",
    ))

    for entity in minimum_entities
        minimum = ParameterValueRef(
            minimum_binding.parameter,
            (parameter_storage_index(plan, minimum_binding.parameter, 1, entity),),
        )
        maximum = ParameterValueRef(
            maximum_binding.parameter,
            (parameter_storage_index(plan, maximum_binding.parameter, 1, entity),),
        )
        push!(constraints, OrderedParameterConstraint(
            minimum, maximum, process_id(named), path, entity,
        ))
    end
    return nothing
end

function _append_quota_factor_constraints!(
    constraints, definition, layout, plan, named, path, factor, refs
)
    factor isa QuotaResponse && _append_quota_constraints!(
        constraints, definition, layout, plan, named, path, refs.slots,
    )
    for (name, subfactor) in pairs(factor_subfactors(factor))
        _append_quota_factor_constraints!(
            constraints,
            definition,
            layout,
            plan,
            named,
            factor_subfactor_path(path, factor, name),
            subfactor,
            getproperty(refs.subfactors, name),
        )
    end
    return nothing
end

function _append_product_constraints!(constraints, definition, named)
    products = process_products(named.process)
    (isnothing(products) || length(products.destinations) == 1) && return nothing
    terms = Tuple(begin
        refs = getproperty(named.binding_refs.products.fractions, product)
        binding = _resolved_slot_bindings(definition, refs).fraction
        ParameterValueRef(binding.parameter, ())
    end for product in keys(products.fractions))
    push!(constraints, ProductFractionConstraint(terms, process_id(named)))
    return nothing
end

"""Build layout-resolved coupled scientific constraints once during model setup."""
function parameter_constraints(
    definition::CanonicalModelDefinition, layout::ModelLayout, plan::ParameterPlan
)
    constraints = Any[]
    for named in values(definition.processes)
        for (name, factor) in pairs(factors(named))
            _append_quota_factor_constraints!(
                constraints,
                definition,
                layout,
                plan,
                named,
                (:factors, name),
                factor,
                getproperty(named.binding_refs.factors, name),
            )
        end
        named.process isa NutrientUptake && _append_quota_constraints!(
            constraints, definition, layout, plan, named, (), named.binding_refs.process,
        )
        _append_product_constraints!(constraints, definition, named)
    end
    return Tuple(constraints)
end

"""Validate formulation-owned scientific domains and coupled constraints on realized values."""
function validate_realized_parameters(
    definition::CanonicalModelDefinition, parameter_values::NamedTuple, constraints::Tuple
)
    _validate_parameter_domains(definition, parameter_values)
    validate_parameter_constraints(parameter_values, constraints)
    return nothing
end
