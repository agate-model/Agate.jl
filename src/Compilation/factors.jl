_factor_input_operand(input::FactorDriver, layout::ModelLayout) =
    input_operand(layout, input.identity)

function _factor_input_operand(input::FactorComponent, layout::ModelLayout)
    target = _scalar_component_target(layout, input.component)
    return input_operand(layout, target)
end

function _factor_element(
    definition::NormalizedModelDefinition,
    named::NamedProcess,
    path::Tuple,
    factor::AbstractFactor,
    layout::ModelLayout,
    plan::ParameterPlan,
    axis_positions::NamedTuple,
)
    input_operands = Tuple(
        _factor_input_operand(input, layout) for input in factor_inputs(factor)
    )
    children = factor_children(factor)
    child_operands = if isempty(children)
        ()
    else
        child_factors = Tuple(
            _factor_element(
                definition,
                named,
                factor_child_path(path, factor, name),
                child,
                layout,
                plan,
                axis_positions,
            )
            for (name, child) in pairs(children)
        )
        (TupleOp(child_factors),)
    end
    slots = parameter_slot_bindings(
        definition,
        named,
        path,
        factor;
        context=factor_parameter_context(factor),
    )
    parameter_operands = Tuple(
        parameter_operand(binding, plan, axis_positions) for binding in values(slots)
    )
    return FactorElement(
        formulation(factor), (input_operands..., child_operands..., parameter_operands...)
    )
end

function _factor_elements(
    definition::NormalizedModelDefinition,
    named::NamedProcess,
    layout::ModelLayout,
    plan::ParameterPlan,
    axis_positions::NamedTuple,
)
    return Tuple(
        _factor_element(
            definition, named, (:factors, name), factor, layout, plan, axis_positions
        )
        for (name, factor) in pairs(factors(named))
    )
end
