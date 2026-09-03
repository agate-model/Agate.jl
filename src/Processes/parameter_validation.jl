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

function _validate_quota_bounds(
    definition::CanonicalModelDefinition,
    parameter_values::NamedTuple,
    named::CanonicalProcess,
    path::Tuple,
    slot_refs::NamedTuple,
)
    slots = _resolved_slot_bindings(definition, slot_refs)
    minimum_binding = slots.minimum_quota
    maximum_binding = slots.maximum_quota
    minimum_values = _parameter_entries(getproperty(parameter_values, minimum_binding.parameter))
    maximum_values = _parameter_entries(getproperty(parameter_values, maximum_binding.parameter))
    for (minimum, maximum) in zip(minimum_values, maximum_values)
        maximum > minimum && continue
        throw(ArgumentError(
            "process :$(process_id(named)) path $path parameter :$(maximum_binding.parameter) " *
            "must be greater than :$(minimum_binding.parameter)=$minimum; got $maximum",
        ))
    end
    return nothing
end

function _validate_quota_factor_bounds(
    definition, parameter_values, named, path, factor, refs
)
    factor isa QuotaResponse && _validate_quota_bounds(
        definition, parameter_values, named, path, refs.slots,
    )
    for (name, subfactor) in pairs(factor_subfactors(factor))
        _validate_quota_factor_bounds(
            definition,
            parameter_values,
            named,
            factor_subfactor_path(path, factor, name),
            subfactor,
            getproperty(refs.subfactors, name),
        )
    end
    return nothing
end

function _validate_quota_bounds(definition, parameter_values)
    for named in values(definition.processes)
        for (name, factor) in pairs(factors(named))
            _validate_quota_factor_bounds(
                definition, parameter_values, named, (:factors, name), factor,
                getproperty(named.binding_refs.factors, name),
            )
        end
        named.process isa NutrientUptake && _validate_quota_bounds(
            definition, parameter_values, named, (), named.binding_refs.process,
        )
    end
    return nothing
end

function _validate_product_fractions(
    definition::CanonicalModelDefinition, parameter_values::NamedTuple
)
    for named in values(definition.processes)
        products = process_products(named.process)
        (isnothing(products) || length(products.destinations) == 1) && continue
        names = keys(products.fractions)
        fractions = Tuple(begin
            refs = getproperty(named.binding_refs.products.fractions, product)
            binding = _resolved_slot_bindings(definition, refs).fraction
            getproperty(parameter_values, binding.parameter)
        end for product in names)

        total = sum(fractions)
        if length(names) == length(products.destinations)
            tolerance = total isa AbstractFloat ? 100 * eps(one(total)) : zero(total)
            isapprox(total, one(total); rtol=zero(total), atol=tolerance) || throw(ArgumentError(
                "explicit product fractions for process :$(process_id(named)) must sum to 1; got $total",
            ))
        else
            remainder = one(total) - total
            zero(remainder) <= remainder <= one(remainder) || throw(ArgumentError(
                "product fractions for process :$(process_id(named)) leave conservative remainder $remainder; expected a value in [0, 1]",
            ))
        end
    end
    return nothing
end

"""Validate formulation-owned scientific domains on fully realized parameter values."""
function validate_realized_parameters(
    definition::CanonicalModelDefinition, parameter_values::NamedTuple
)
    _validate_parameter_domains(definition, parameter_values)
    _validate_quota_bounds(definition, parameter_values)
    _validate_product_fractions(definition, parameter_values)
    return nothing
end
