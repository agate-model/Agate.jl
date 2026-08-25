function process_fluxes(
    named::NamedProcess{P},
    definition::NormalizedModelDefinition,
    layout::ComponentLayout,
    context::CommunityContext,
) where {P<:Remineralization}
    process = named.process
    destination = _scalar_component_target(layout, process.destination)
    fluxes = Any[]

    for source_component in process.sources
        source = _scalar_component_target(layout, source_component)
        rate_binding = parameter_slot_bindings(
            definition,
            named,
            (),
            process;
            context=(source=source_component,),
        ).rate
        rate = RateElement(
            process.formulation,
            (TracerOp{source}(), parameter_operand(rate_binding)),
        )
        push!(
            fluxes,
            FluxSpec(process_id(named), source, rate, Weight{-1}()),
            FluxSpec(process_id(named), destination, rate, Weight{1}()),
        )
    end
    return Tuple(fluxes)
end
