function _product_fraction_binding(
    context::CompileContext,
    named::NamedProcess,
    products::Products,
    product::Symbol,
)
    return parameter_slot_bindings(
        context.definition,
        named,
        product_path(named.process),
        products;
        context=(product=product,),
    ).fraction
end

function _product_fraction_operand(
    context::CompileContext,
    named::NamedProcess,
    products::Products,
    product::Symbol,
)
    if product !== products.balanced
        hasproperty(products.fractions, product) || return nothing
        return parameter_operand(
            _product_fraction_binding(context, named, products, product), context.plan
        )
    end

    explicit_products = Tuple(
        name for name in keys(products.targets)
        if name !== products.balanced && hasproperty(products.fractions, name)
    )
    isempty(explicit_products) && return nothing
    operands = Tuple(
        parameter_operand(
            _product_fraction_binding(context, named, products, name), context.plan
        )
        for name in explicit_products
    )
    return RemainderOp(operands)
end

function _product_ratio_binding(
    context::CompileContext,
    named::NamedProcess,
    products::Products,
    currency::Symbol,
)
    stoichiometry = products.stoichiometry
    isnothing(stoichiometry) && return nothing
    currency === stoichiometry.reference && return nothing
    return parameter_slot_bindings(
        context.definition,
        named,
        (product_path(named.process)..., :stoichiometry),
        stoichiometry;
        context=(currency=currency,),
    ).ratio
end

_product_operand_tuple(::Nothing) = ()
_product_operand_tuple(operand) = (operand,)

function _product_weight(
    fraction,
    ratio::Union{Nothing,ParameterBinding},
    plan::ParameterPlan;
    suffix::Tuple=(),
)
    ratio_operand = isnothing(ratio) ? nothing : parameter_operand(ratio, plan)
    operands = (
        _product_operand_tuple(fraction)...,
        _product_operand_tuple(ratio_operand)...,
        suffix...,
    )
    return Weight{1}(operands)
end

function _product_fluxes(
    named::NamedProcess,
    product_targets::NamedTuple,
    context::CompileContext,
    rate::RateElement;
    suffix::Tuple=(),
)
    products = process_products(named.process)
    fluxes = Any[]
    for (product, targets) in pairs(product_targets)
        fraction = _product_fraction_operand(context, named, products, product)
        for (currency, component) in pairs(targets)
            target = _scalar_component_target(context.layout, component)
            ratio = _product_ratio_binding(context, named, products, currency)
            push!(
                fluxes,
                FluxSpec(
                    target, rate,
                    _product_weight(fraction, ratio, context.plan; suffix),
                ),
            )
        end
    end
    return Tuple(fluxes)
end
