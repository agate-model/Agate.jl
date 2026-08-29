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

_factor_process_operands(
    ::AbstractFactor, ::CompileContext, ::NamedProcess, ::NamedTuple
) = ()

_factor_inputs(factor::AbstractFactor, ::NamedProcess) = factor_inputs(factor)
function _factor_inputs(factor::QuotaResponse, named::NamedProcess)
    reference = only(named.facts.population_states)
    variable = PopulationStateRef(reference.population, factor.variable_state)
    return (FactorPopulationState(variable), FactorPopulationState(reference))
end

function _factor_process_operands(
    ::Light,
    context::CompileContext,
    named::NamedProcess,
    axis_positions::NamedTuple,
)
    ref = named.binding_refs.process.maximum_rate
    return (parameter_operand(ref, context, axis_positions),)
end

function _factor_element(
    context::CompileContext,
    named::NamedProcess,
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
    children = factor_children(factor)
    child_operands = if isempty(children)
        ()
    else
        child_factors = Tuple(
            _factor_element(
                context,
                named,
                child,
                getproperty(refs.children, name),
                axis_positions,
            )
            for (name, child) in pairs(children)
        )
        (TupleOp(child_factors),)
    end
    parameter_operands = Tuple(
        parameter_operand(ref, context, axis_positions) for ref in values(refs.slots)
    )
    return FactorElement(
        formulation(factor),
        (input_operands..., process_operands..., child_operands..., parameter_operands...),
    )
end

function _factor_elements(
    context::CompileContext,
    named::NamedProcess,
    axis_positions::NamedTuple,
)
    return Tuple(
        _factor_element(
            context, named, factor, getproperty(named.binding_refs.factors, name), axis_positions
        )
        for (name, factor) in pairs(factors(named))
    )
end
