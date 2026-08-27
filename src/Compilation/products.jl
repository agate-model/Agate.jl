function _product_fraction_ref(named::NamedProcess, product::Symbol)
    refs = named.binding_refs.products.fractions
    hasproperty(refs, product) || return nothing
    return getproperty(refs, product).fraction
end
function _product_fraction_operand(
    context::CompileContext,
    named::NamedProcess,
    products::Products,
    product::Symbol,
)
    if product !== products.balanced
        ref = _product_fraction_ref(named, product)
        isnothing(ref) && return nothing
        return parameter_operand(ref, context)
    end

    fraction_refs = named.binding_refs.products.fractions
    explicit_products = Tuple(
        name for name in keys(products.targets)
        if name !== products.balanced && hasproperty(fraction_refs, name)
    )
    isempty(explicit_products) && return nothing
    operands = Tuple(
        parameter_operand(_product_fraction_ref(named, name), context)
        for name in explicit_products
    )
    return RemainderOp(operands)
end

function _product_ratio_ref(named::NamedProcess, products::Products, currency::Symbol)
    stoichiometry = products.stoichiometry
    isnothing(stoichiometry) && return nothing
    currency === stoichiometry.reference && return nothing
    return getproperty(named.binding_refs.products.stoichiometry, currency).ratio
end
_product_operand_tuple(::Nothing) = ()
_product_operand_tuple(operand) = (operand,)

function _product_weight(
    fraction,
    ratio::Union{Nothing,Int},
    context::CompileContext;
    suffix::Tuple=(),
)
    ratio_operand = isnothing(ratio) ? nothing : parameter_operand(ratio, context)
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
            ratio = _product_ratio_ref(named, products, currency)
            push!(
                fluxes,
                FluxSpec(
                    target, rate,
                    _product_weight(fraction, ratio, context; suffix),
                ),
            )
        end
    end
    return Tuple(fluxes)
end
