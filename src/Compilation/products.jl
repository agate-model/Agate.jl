function _product_fraction_ref(named::CanonicalProcess, product::Symbol)
    refs = named.binding_refs.products.fractions
    hasproperty(refs, product) || return nothing
    return getproperty(refs, product).fraction
end

function _product_fraction_operand(
    context::CompileContext,
    named::CanonicalProcess,
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
    return ComplementOp(operands)
end

"""Return the implicit material ratio operand for this plankton entity, if needed."""
function _material_ratio_operand(context::CompileContext, plankton::Symbol, element::Symbol, position)
    ratios = getproperty(context.definition.composition_refs, plankton)
    hasproperty(ratios, element) || return ()
    ref = getproperty(ratios, element)
    return (parameter_operand(ref, context, (plankton=position,)),)
end

_product_operand_tuple(::Nothing) = ()
_product_operand_tuple(operand) = (operand,)

function _product_weight(fraction; suffix::Tuple=())
    operands = (_product_operand_tuple(fraction)..., suffix...)
    return Weight{1}(operands)
end

function _product_fluxes_for_element(
    named::CanonicalProcess,
    product_targets::NamedTuple,
    context::CompileContext,
    rate::RateOp,
    element::Symbol;
    suffix::Tuple=(),
)
    fluxes = Any[]
    for (product, targets) in pairs(product_targets)
        fraction = _product_fraction_operand(context, named, product)
        target = _scalar_component_target(context.layout, getproperty(targets, element))
        push!(
            fluxes,
            FluxSpec(target, rate, _product_weight(fraction; suffix)),
        )
    end
    return Tuple(fluxes)
end
