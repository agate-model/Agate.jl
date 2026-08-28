function process_fluxes(
    named::NamedProcess{P}, context::CompileContext
) where {P<:Remineralization}
    process = named.process
    layout = context.layout
    destination = _scalar_component_target(layout, process.destination)
    fluxes = Any[]

    for source_component in process.sources
        slots = getproperty(named.binding_refs.process, source_component)
        source = _scalar_component_target(layout, source_component)
        rate_ref = slots.rate
        rate = RateElement(
            process.formulation,
            (input_operand(layout, source), parameter_operand(rate_ref, context)),
        )
        push!(
            fluxes,
            FluxSpec(source, rate, Weight{-1}()),
            FluxSpec(destination, rate, Weight{1}()),
        )
    end
    return Tuple(fluxes)
end
