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

function _routing_fluxes(
    named::NamedProcess,
    definition::NormalizedModelDefinition,
    routing::ProductRouting{DirectRouting},
    layout::ComponentLayout,
    rate::RateElement;
    suffix::Tuple=(),
)
    target = _scalar_component_target(layout, routing.retained)
    return (FluxSpec(process_id(named), target, rate, Weight{1}(suffix)),)
end

function _routing_fluxes(
    named::NamedProcess,
    definition::NormalizedModelDefinition,
    routing::ProductRouting{PartitionRouting},
    layout::ComponentLayout,
    rate::RateElement;
    suffix::Tuple=(),
)
    fraction = _routing_fraction_binding(definition, named, routing)
    retained = FluxSpec(
        process_id(named),
        _scalar_component_target(layout, routing.retained),
        rate,
        _routing_weight(fraction, :retained; suffix),
    )
    exported = FluxSpec(
        process_id(named),
        _scalar_component_target(layout, routing.exported),
        rate,
        _routing_weight(fraction, :exported; suffix),
    )
    return (retained, exported)
end

function _routing_fluxes(
    named::NamedProcess,
    definition::NormalizedModelDefinition,
    routing::ProductRouting{DOMPOMRouting},
    layout::ComponentLayout,
    rate::RateElement;
    suffix::Tuple=(),
)
    fraction = _routing_fraction_binding(definition, named, routing)
    fluxes = ()
    for route in (:DOM, :POM)
        pool = getproperty(routing.pools, route)
        for currency in keys(pool)
            target = _scalar_component_target(layout, getproperty(pool, currency))
            ratio = _routing_ratio_binding(definition, named, routing, currency)
            weight = _routing_weight(fraction, route; ratio, suffix)
            fluxes = (fluxes..., FluxSpec(process_id(named), target, rate, weight))
        end
    end
    return fluxes
end

function _routing_fluxes(
    named::NamedProcess,
    definition::NormalizedModelDefinition,
    routing::ProductRouting,
    layout::ComponentLayout,
    rate::RateElement;
    suffix::Tuple=(),
)
    throw(ArgumentError("unsupported routing formulation $(typeof(routing.formulation))"))
end
