"""Validated process instance in canonical model state.

`semantic_facts` contains setup-only scientific decisions that compilation may trust. `binding_refs`
contains dense references into the canonical model's ordered parameter-binding tuple, arranged
alongside the process/factor/product structure so lowering never reconstructs scientific paths.
"""
struct CanonicalProcess{Process<:AbstractProcess,SemanticFacts,BindingRefs}
    id::Symbol
    process::Process
    semantic_facts::SemanticFacts
    binding_refs::BindingRefs
end

CanonicalProcess(id::Symbol, process::Process) where {Process<:AbstractProcess} =
    CanonicalProcess(id, process, NamedTuple(), NamedTuple())
CanonicalProcess(
    id::Symbol, process::Process, semantic_facts::SemanticFacts
) where {Process<:AbstractProcess,SemanticFacts} =
    CanonicalProcess(id, process, semantic_facts, NamedTuple())

"""Return the stable identity of a canonical process."""
process_id(process::CanonicalProcess) = process.id
formulation(process::CanonicalProcess) = formulation(process.process)
factors(process::CanonicalProcess) = factors(process.process)
participants(process::CanonicalProcess) = participants(process.process)

# Dense parameter-binding reference collection

function _resolve_slot_qualifier(slot::ParameterSlot, context::NamedTuple)
    isnothing(slot.qualify) && return nothing
    hasproperty(context, slot.qualify) || return nothing
    value = getproperty(context, slot.qualify)
    value isa Symbol || throw(
        ArgumentError("parameter slot qualifier :$(slot.qualify) must identify one Symbol"),
    )
    return Qualifier(slot.qualify, value)
end

function _binding_axis_components(
    named::CanonicalProcess, axes::Tuple, qualifier
)
    process_participants = participants(named)
    return map(axes) do axis
        qualifier isa Qualifier && qualifier.axis === axis && return (qualifier.value,)
        hasproperty(process_participants, axis) || throw(ArgumentError(
            "parameter applicability axis :$axis is not a participant role of process :$(process_id(named))",
        ))
        getproperty(process_participants, axis)
    end
end

function _parameter_slot_metadata(
    named::CanonicalProcess,
    path::Tuple,
    slot::ParameterSlot,
    qualifier,
)
    return (;
        process=process_id(named),
        path,
        slot=slot.name,
        qualifier,
        axes=slot.axes,
        domain=slot.domain,
        axis_components=_binding_axis_components(named, slot.axes, qualifier),
    )
end

function _emit_parameter_slots!(
    uses::Vector{Any},
    seen::Set{Any},
    named::CanonicalProcess,
    path::Tuple,
    node;
    context::NamedTuple=NamedTuple(),
)
    slots = parameter_slots(_parameter_slot_source(node))
    names = Tuple(slot.name for slot in slots)
    bindings = _resolve_authored_bindings(node)
    refs = Vector{Int}(undef, length(slots))
    for i in eachindex(slots)
        slot = slots[i]
        qualifier = _resolve_slot_qualifier(slot, context)
        metadata = _parameter_slot_metadata(named, path, slot, qualifier)
        qualifier_key = isnothing(qualifier) ? nothing : (qualifier.axis, qualifier.value)
        identity = (metadata.process, metadata.path, metadata.slot, qualifier_key)
        identity in seen && throw(
            ArgumentError("canonical processes declare duplicate parameter binding key $identity"),
        )
        push!(seen, identity)
        parameter, explicit, binding_qualifiers = _resolve_binding_value(bindings, slot, qualifier)
        push!(uses, (; metadata..., parameter, explicit, binding_qualifiers))
        refs[i] = length(uses)
    end
    return NamedTuple{names}(Tuple(refs))
end

function _visit_factor_slots!(
    uses::Vector{Any}, seen::Set{Any}, named::CanonicalProcess, path::Tuple, factor::AbstractFactor
)
    slots = _emit_parameter_slots!(uses, seen, named, path, factor)
    subfactors = _resolve_factor_subfactors(factor)
    names = keys(subfactors)
    subfactor_refs = NamedTuple{names}(Tuple(
        _visit_factor_slots!(
            uses, seen, named, factor_subfactor_path(path, factor, name), subfactor
        )
        for (name, subfactor) in pairs(subfactors)
    ))
    return (; slots, subfactors=subfactor_refs)
end

function _visit_product_slots!(
    uses::Vector{Any}, seen::Set{Any}, named::CanonicalProcess, path::Tuple, products::Products
)
    fraction_names = keys(products.fractions)
    fractions = NamedTuple{fraction_names}(Tuple(
        _emit_parameter_slots!(
            uses,
            seen,
            named,
            path,
            products;
            context=(product=product,),
        )
        for product in fraction_names
    ))

    stoichiometry = products.stoichiometry
    isnothing(stoichiometry) && return (; fractions, stoichiometry=NamedTuple())
    elements = Tuple(
        element for element in keys(first(values(products.destinations)))
        if element !== stoichiometry.reference_element
    )
    ratios = NamedTuple{elements}(Tuple(
        _emit_parameter_slots!(
            uses,
            seen,
            named,
            (path..., :stoichiometry),
            stoichiometry;
            context=(element=element,),
        )
        for element in elements
    ))
    return (; fractions, stoichiometry=ratios)
end

function _visit_process_slots!(uses::Vector{Any}, seen::Set{Any}, named::CanonicalProcess)
    process = named.process
    process_refs = if process isa Mortality
        names = process.plankton
        NamedTuple{names}(Tuple(
            _emit_parameter_slots!(
                uses, seen, named, (), process; context=(plankton=plankton,)
            )
            for plankton in names
        ))
    elseif process isa Remineralization
        names = process.sources
        NamedTuple{names}(Tuple(
            _emit_parameter_slots!(
                uses, seen, named, (), process; context=(source=source,)
            )
            for source in names
        ))
    else
        _emit_parameter_slots!(uses, seen, named, (), process)
    end

    process_factors = factors(process)
    factor_names = keys(process_factors)
    factor_refs = NamedTuple{factor_names}(Tuple(
        _visit_factor_slots!(uses, seen, named, (:factors, name), factor)
        for (name, factor) in pairs(process_factors)
    ))

    products = process_products(process)
    product_refs = isnothing(products) ? nothing :
        _visit_product_slots!(uses, seen, named, product_path(process), products)

    stoichiometry_refs = NamedTuple()
    if process isa Growth && !isnothing(process.stoichiometry)
        elements = keys(process.additional_resources)
        stoichiometry_refs = NamedTuple{elements}(Tuple(
            _emit_parameter_slots!(
                uses,
                seen,
                named,
                (:stoichiometry,),
                process.stoichiometry;
                context=(element=element,),
            )
            for element in elements
        ))
    end
    return (;
        process=process_refs,
        factors=factor_refs,
        products=product_refs,
        stoichiometry=stoichiometry_refs,
    )
end
