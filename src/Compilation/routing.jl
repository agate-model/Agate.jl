function _realize_routing_targets(pool::NamedTuple, layout::ComponentLayout)
    names = keys(pool)
    values = Tuple(_scalar_component_target(layout, getproperty(pool, name)) for name in names)
    return NamedTuple{names}(values)
end

function _realize_dom_pom_routing(routing::ProductRouting, layout::ComponentLayout)
    routing.formulation isa DOMPOMRouting || throw(
        ArgumentError("routing is not a DOM/POM routing formulation"),
    )
    return (
        DOM=_realize_routing_targets(routing.pools.DOM, layout),
        POM=_realize_routing_targets(routing.pools.POM, layout),
    )
end

function _routing_fraction_binding(
    definition::NormalizedModelDefinition, named::NamedProcess, routing::ProductRouting
)
    slots = parameter_slot_bindings(definition, named, (:routing,), routing.formulation)
    routing.formulation isa PartitionRouting && return slots.export_fraction
    routing.formulation isa DOMPOMRouting && return slots.POM_fraction
    throw(ArgumentError("unsupported routing formulation $(typeof(routing.formulation))"))
end

function _routing_ratio_binding(
    definition::NormalizedModelDefinition,
    named::NamedProcess,
    routing::ProductRouting,
    currency::Symbol,
)
    currency === routing.stoichiometry.reference && return nothing
    return parameter_slot_bindings(
        definition,
        named,
        (:routing, :stoichiometry),
        routing.stoichiometry;
        context=(currency=currency,),
    ).ratio
end

function _fraction_operand(binding::ParameterBinding, route::Symbol)
    operand = parameter_operand(binding)
    route in (:retained, :DOM) && return ComplementOp(operand)
    route in (:exported, :POM) && return operand
    throw(ArgumentError("unsupported routing partition :$route"))
end

_ratio_operands(::Nothing) = ()
_ratio_operands(binding::ParameterBinding) = (parameter_operand(binding),)

function _routing_weight(
    fraction::ParameterBinding,
    route::Symbol;
    ratio::Union{Nothing,ParameterBinding}=nothing,
    suffix::Tuple=(),
)
    operands = (_fraction_operand(fraction, route), _ratio_operands(ratio)..., suffix...)
    return Weight{1}(operands)
end
