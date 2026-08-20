_factor_input_operand(input::FactorDriver, layout::ComponentLayout) =
    TracerOp{input.identity}()

function _factor_input_operand(input::FactorComponent, layout::ComponentLayout)
    target = _scalar_component_target(layout, input.component)
    return TracerOp{target}()
end

function _factor_parameter_operand(binding::ParameterBinding, axis_indices::NamedTuple)
    indices = Tuple(
        begin
            hasproperty(axis_indices, axis) || throw(ArgumentError(
                "factor parameter axis :$axis has no realized runtime index",
            ))
            getproperty(axis_indices, axis)
        end
        for axis in binding.requirement.axes
    )
    return parameter_operand(binding, indices...)
end

function _factor_parameter_operands(bindings::NamedTuple, axis_indices::NamedTuple)
    return Tuple(
        _factor_parameter_operand(binding, axis_indices) for binding in values(bindings)
    )
end

function _factor_element(
    definition::NormalizedModelDefinition,
    named::NamedProcess,
    path::Tuple,
    factor::AbstractFactor,
    layout::ComponentLayout,
    axis_indices::NamedTuple,
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
                axis_indices,
            )
            for (name, child) in pairs(children)
        )
        (TupleOp(child_factors),)
    end
    slots = parameter_slot_bindings(
        definition,
        named,
        path,
        formulation(factor);
        context=factor_parameter_context(factor),
    )
    parameter_operands = _factor_parameter_operands(slots, axis_indices)
    return FactorElement(
        formulation(factor), (input_operands..., child_operands..., parameter_operands...)
    )
end

function _factor_elements(
    definition::NormalizedModelDefinition,
    named::NamedProcess,
    layout::ComponentLayout,
    axis_indices::NamedTuple,
)
    return Tuple(
        _factor_element(
            definition, named, (:factors, name), factor, layout, axis_indices
        )
        for (name, factor) in pairs(factors(named))
    )
end
