_factor_input_operand(input::FactorDriver, layout::ModelLayout, axis_positions::NamedTuple) =
    input_operand(layout, input.identity)

function _factor_input_operand(
    input::FactorComponent, layout::ModelLayout, axis_positions::NamedTuple
)
    target = _scalar_component_target(layout, input.component)
    return input_operand(layout, target)
end

function _factor_input_operand(
    input::FactorPlanktonState, layout::ModelLayout, axis_positions::NamedTuple
)
    hasproperty(axis_positions, :plankton) || throw(ArgumentError(
        "plankton-state factor input requires a realized :plankton axis position",
    ))
    position = axis_positions.plankton
    position.component === input.reference.plankton || throw(ArgumentError(
        "plankton-state factor input :$(input.reference.plankton).$(input.reference.state) " *
        "must reference the current logical plankton :$(position.component)",
    ))
    tracer = state_tracer(layout, input.reference, position.component_index)
    return input_operand(layout, tracer)
end

_factor_process_operands(
    ::AbstractFactor, ::CompileContext, ::CanonicalProcess, ::NamedTuple
) = ()

_factor_inputs(factor::AbstractFactor, ::CanonicalProcess) = factor_inputs(factor)
function _factor_inputs(factor::QuotaResponse, named::CanonicalProcess)
    reference = only(named.semantic_facts.plankton_states)
    variable = PlanktonStateRef(reference.plankton, factor.variable_state)
    return (FactorPlanktonState(variable), FactorPlanktonState(reference))
end

function _factor_process_operands(
    ::Light,
    context::CompileContext,
    named::CanonicalProcess,
    axis_positions::NamedTuple,
)
    ref = named.binding_refs.process.maximum_rate
    return (parameter_operand(ref, context, axis_positions),)
end

function _factor_op(
    context::CompileContext,
    named::CanonicalProcess,
    factor::AbstractFactor,
    refs,
    axis_positions::NamedTuple,
)
    layout = context.layout
    input_operands = Tuple(
        _factor_input_operand(input, layout, axis_positions) for input in _factor_inputs(factor, named)
    )
    process_operands = _factor_process_operands(
        factor, context, named, axis_positions
    )
    subfactors = factor_subfactors(factor)
    subfactor_operands = if isempty(subfactors)
        ()
    else
        subfactor_ops = Tuple(
            _factor_op(
                context,
                named,
                subfactor,
                getproperty(refs.subfactors, name),
                axis_positions,
            )
            for (name, subfactor) in pairs(subfactors)
        )
        (TupleOp(subfactor_ops),)
    end
    parameter_operands = Tuple(
        parameter_operand(ref, context, axis_positions) for ref in values(refs.slots)
    )
    return FactorOp(
        formulation(factor),
        (input_operands..., process_operands..., subfactor_operands..., parameter_operands...),
    )
end

function _factor_ops(
    context::CompileContext,
    named::CanonicalProcess,
    axis_positions::NamedTuple,
)
    return Tuple(
        _factor_op(
            context, named, factor, getproperty(named.binding_refs.factors, name), axis_positions
        )
        for (name, factor) in pairs(factors(named))
    )
end
