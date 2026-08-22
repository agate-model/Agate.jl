_factor_input_operand(input::FactorDriver, layout::ComponentLayout) =
    TracerOp{input.identity}()

function _factor_input_operand(input::FactorComponent, layout::ComponentLayout)
    target = _scalar_component_target(layout, input.component)
    return TracerOp{target}()
end

function _factor_parameter_operand(
    binding::ParameterBinding, context::CommunityContext, axis_positions::NamedTuple
)
    return parameter_operand(binding, context, axis_positions)
end

function _factor_parameter_operands(
    bindings::NamedTuple, context::CommunityContext, axis_positions::NamedTuple
)
    return Tuple(
        _factor_parameter_operand(binding, context, axis_positions)
        for binding in values(bindings)
    )
end

function _factor_element(
    definition::NormalizedModelDefinition,
    named::NamedProcess,
    path::Tuple,
    factor::AbstractFactor,
    layout::ComponentLayout,
    context::CommunityContext,
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
                context,
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
    parameter_operands = _factor_parameter_operands(slots, context, axis_positions)
    return FactorElement(
        formulation(factor), (input_operands..., child_operands..., parameter_operands...)
    )
end

function _factor_elements(
    definition::NormalizedModelDefinition,
    named::NamedProcess,
    layout::ComponentLayout,
    context::CommunityContext,
    axis_positions::NamedTuple,
)
    return Tuple(
        _factor_element(
            definition, named, (:factors, name), factor, layout, context, axis_positions
        )
        for (name, factor) in pairs(factors(named))
    )
end
