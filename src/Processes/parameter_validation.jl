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
        fractions = Tuple(begin
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
        end for product in names)

        total = sum(fractions)
        if length(names) == length(products.targets)
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

"""Validate numeric scientific constraints on fully realized parameter values."""
function validate_realized_parameters(
    definition::CanonicalModelDefinition,
    layout::ModelLayout,
    plan::ParameterPlan,
    parameter_values::NamedTuple,
)
    _validate_product_fractions(definition, parameter_values)
    _validate_quota_science(definition, layout, plan, parameter_values)
    return nothing
end
