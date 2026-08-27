_factor_input_operand(input::FactorDriver, layout::ModelLayout, axis_positions::NamedTuple) =
    input_operand(layout, input.identity)

function _factor_input_operand(
    input::FactorComponent, layout::ModelLayout, axis_positions::NamedTuple
)
    target = _scalar_component_target(layout, input.component)
    return input_operand(layout, target)
end

function _factor_input_operand(
    input::FactorPopulationState, layout::ModelLayout, axis_positions::NamedTuple
)
    hasproperty(axis_positions, :population) || throw(ArgumentError(
        "population-state factor input requires a realized :population axis position",
    ))
    local_index = axis_positions.population.local_index
    tracer = state_tracer(layout, input.reference, local_index)
    return input_operand(layout, tracer)
end

function _factor_element(
    context::CompileContext,
    named::NamedProcess,
    path::Tuple,
    factor::AbstractFactor,
    axis_positions::NamedTuple,
)
    layout = context.layout
    input_operands = Tuple(
        _factor_input_operand(input, layout, axis_positions) for input in factor_inputs(factor)
    )
    children = factor_children(factor)
    child_operands = if isempty(children)
        ()
    else
        child_factors = Tuple(
            _factor_element(
                context,
                named,
                factor_child_path(path, factor, name),
                child,
                axis_positions,
            )
            for (name, child) in pairs(children)
        )
        (TupleOp(child_factors),)
    end
    slots = parameter_slot_bindings(
        context.definition,
        named,
        path,
        factor;
        context=factor_parameter_context(factor),
    )
    parameter_operands = Tuple(
        parameter_operand(binding, context.plan, axis_positions) for binding in values(slots)
    )
    return FactorElement(
        formulation(factor), (input_operands..., child_operands..., parameter_operands...)
    )
end

function _factor_elements(
    context::CompileContext,
    named::NamedProcess,
    axis_positions::NamedTuple,
)
    return Tuple(
        _factor_element(context, named, (:factors, name), factor, axis_positions)
        for (name, factor) in pairs(factors(named))
    )
end
