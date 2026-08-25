function _product_fraction_binding(
    definition::NormalizedModelDefinition,
    named::NamedProcess,
    products::Products,
    route::Symbol,
)
    return parameter_slot_bindings(
        definition, named, product_path(named.process), products; context=(product=route,)
    ).fraction
end

function _product_fraction_operand(
    definition::NormalizedModelDefinition,
    named::NamedProcess,
    products::Products,
    route::Symbol,
)
    hasproperty(products.fractions, route) && return parameter_operand(
        _product_fraction_binding(definition, named, products, route)
    )
    route === products.inferred || return nothing

    explicit_routes = Tuple(keys(products.fractions))
    isempty(explicit_routes) && return nothing
    operands = Tuple(
        parameter_operand(_product_fraction_binding(definition, named, products, name))
        for name in explicit_routes
    )
    length(operands) == 1 && return ComplementOp(only(operands))
    return OneMinusSumOp(operands)
end

function _product_fluxes(
    named::NamedProcess,
    definition::NormalizedModelDefinition,
    products::Products,
    layout::ComponentLayout,
    rate::RateElement;
    suffix::Tuple=(),
)
    fluxes = Any[]
    for route in keys(products.targets)
        target = _scalar_component_target(layout, getproperty(products.targets, route))
        fraction = _product_fraction_operand(definition, named, products, route)
        operands = isnothing(fraction) ? suffix : (fraction, suffix...)
        push!(fluxes, FluxSpec(process_id(named), target, rate, Weight{1}(operands)))
    end
    return Tuple(fluxes)
end
