"""Resolved scalar parameters for one DOM/POM product-routing rule."""
struct DOMPOMRoutingBinding{R}
    fraction::Symbol
    ratios::R
end

function _dom_pom_routing_binding(
    definition::NormalizedModelDefinition, named::NamedProcess, routing::ProductRouting
)
    routing.formulation isa DOMPOMRouting || throw(
        ArgumentError("routing is not a DOM/POM routing formulation"),
    )
    fraction = parameter_name(
        definition, _parameter_requirement(named, (:routing,), :POM_fraction)
    )
    reference = routing.stoichiometry.reference
    currencies = Tuple(currency for currency in keys(routing.pools.DOM) if currency !== reference)
    ratio_values = Tuple(
        parameter_name(
            definition,
            _parameter_requirement(
                named,
                (:routing, :stoichiometry),
                :ratio;
                qualifier=(currency=currency,),
            ),
        ) for currency in currencies
    )
    return DOMPOMRoutingBinding(fraction, NamedTuple{currencies}(ratio_values))
end

"""Concrete tracer targets for DOM/POM routing, keyed by product currency."""
struct DOMPOMRoutingTopology{D,P}
    DOM::D
    POM::P
    reference::Symbol
end

function _realize_routing_targets(pool::NamedTuple, layout::ComponentLayout)
    names = keys(pool)
    values = Tuple(_scalar_component_target(layout, getproperty(pool, name)) for name in names)
    return NamedTuple{names}(values)
end

function _realize_dom_pom_routing(routing::ProductRouting, layout::ComponentLayout)
    routing.formulation isa DOMPOMRouting || throw(
        ArgumentError("routing is not a DOM/POM routing formulation"),
    )
    return DOMPOMRoutingTopology(
        _realize_routing_targets(routing.pools.DOM, layout),
        _realize_routing_targets(routing.pools.POM, layout),
        routing.stoichiometry.reference,
    )
end

function _routing_ratio_parameter(
    routing::DOMPOMRoutingTopology,
    binding::DOMPOMRoutingBinding,
    currency::Symbol,
)
    currency === routing.reference && return nothing
    return getproperty(binding.ratios, currency)
end

function _fraction_operand(fraction::Symbol, route::Symbol)
    operand = ScalarParamOp{fraction}()
    route in (:retained, :DOM) && return ComplementOp(operand)
    route in (:exported, :POM) && return operand
    throw(ArgumentError("unsupported routing partition :$route"))
end

_ratio_operands(::Nothing) = ()
_ratio_operands(ratio::Symbol) = (ScalarParamOp{ratio}(),)

function _routing_weight(
    fraction::Symbol,
    route::Symbol;
    ratio::Union{Nothing,Symbol}=nothing,
    suffix::Tuple=(),
)
    operands = (_fraction_operand(fraction, route), _ratio_operands(ratio)..., suffix...)
    return Weight{1}(operands)
end
