function _routing_fraction_binding(
    definition::NormalizedModelDefinition, named::NamedProcess, routing::ProductRouting
)
    return parameter_slot_bindings(definition, named, (:routing,), routing).POM_fraction
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

function _routing_fraction_operand(binding::ParameterBinding, route::Symbol)
    operand = parameter_operand(binding)
    route === :DOM && return ComplementOp(operand)
    route === :POM && return operand
    throw(ArgumentError("unsupported DOM/POM route :$route"))
end

_ratio_operands(::Nothing) = ()
_ratio_operands(binding::ParameterBinding) = (parameter_operand(binding),)

function _routing_weight(
    fraction::ParameterBinding,
    route::Symbol;
    ratio::Union{Nothing,ParameterBinding}=nothing,
    suffix::Tuple=(),
)
    operands = (_routing_fraction_operand(fraction, route), _ratio_operands(ratio)..., suffix...)
    return Weight{1}(operands)
end

function _routing_fluxes(
    named::NamedProcess,
    definition::NormalizedModelDefinition,
    routing::ProductRouting,
    layout::ComponentLayout,
    rate::RateElement;
    suffix::Tuple=(),
)
    fraction = _routing_fraction_binding(definition, named, routing)
    fluxes = Any[]
    for route in (:DOM, :POM)
        pool = getproperty(routing.pools, route)
        for currency in keys(pool)
            target = _scalar_component_target(layout, getproperty(pool, currency))
            ratio = _routing_ratio_binding(definition, named, routing, currency)
            weight = _routing_weight(fraction, route; ratio, suffix)
            push!(fluxes, FluxSpec(process_id(named), target, rate, weight))
        end
    end
    return Tuple(fluxes)
end
