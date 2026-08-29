function _product_fraction_ref(named::NamedProcess, product::Symbol)
    refs = named.binding_refs.products.fractions
    hasproperty(refs, product) || return nothing
    return getproperty(refs, product).fraction
end

function _product_fraction_operand(
    context::CompileContext,
    named::NamedProcess,
    product::Symbol,
)
    ref = _product_fraction_ref(named, product)
    isnothing(ref) || return parameter_operand(ref, context)

    fraction_refs = named.binding_refs.products.fractions
    isempty(fraction_refs) && return nothing
    operands = Tuple(
        parameter_operand(getproperty(fraction_refs, name).fraction, context)
        for name in keys(fraction_refs)
    )
    return RemainderOp(operands)
end

function _product_ratio_ref(named::NamedProcess, products::Products, element::Symbol)
    stoichiometry = products.stoichiometry
    isnothing(stoichiometry) && return nothing
    element === stoichiometry.reference_element && return nothing
    return getproperty(named.binding_refs.products.stoichiometry, element).ratio
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
    rate::RateOp;
    suffix::Tuple=(),
)
    products = process_products(named.process)
    fluxes = Any[]
    for (product, targets) in pairs(product_targets)
        fraction = _product_fraction_operand(context, named, product)
        for (element, component) in pairs(targets)
            target = _scalar_component_target(context.layout, component)
            ratio = _product_ratio_ref(named, products, element)
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


function _product_fluxes_for_element(
    named::NamedProcess,
    product_targets::NamedTuple,
    context::CompileContext,
    rate::RateOp,
    element::Symbol;
    suffix::Tuple=(),
)
    products = process_products(named.process)
    fluxes = Any[]
    for (product, targets) in pairs(product_targets)
        fraction = _product_fraction_operand(context, named, product)
        target = _scalar_component_target(context.layout, getproperty(targets, element))
        push!(
            fluxes,
            FluxSpec(target, rate, _product_weight(fraction, nothing, context; suffix)),
        )
    end
    return Tuple(fluxes)
end
