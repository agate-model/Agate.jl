"""Realized linear remineralization topology over scalar source and destination tracers."""
struct RemineralizationTopology{S,D}
    source_tracers::S
    destination_target::D
end

function realize_process_topology(
    named::NamedProcess{P}, layout::ComponentLayout, context::CommunityContext
) where {P<:Remineralization}
    process = named.process
    process.formulation isa LinearRemineralization || throw(
        ArgumentError("unsupported remineralization formulation $(typeof(process.formulation))"),
    )
    length(process.destinations) == 1 || throw(
        ArgumentError("linear remineralization currently requires one destination"),
    )
    sources = Tuple(_scalar_component_target(layout, source) for source in process.sources)
    destination = _scalar_component_target(layout, only(process.destinations))
    return RemineralizationTopology(sources, destination)
end

function process_fluxes(
    named::NamedProcess{P},
    topology::RemineralizationTopology,
    definition::NormalizedModelDefinition,
) where {P<:Remineralization}
    process = named.process
    fluxes = ()
    for i in eachindex(topology.source_tracers)
        source = topology.source_tracers[i]
        component = process.sources[i]
        rate_binding = parameter_slot_bindings(
            definition,
            named,
            (),
            process.formulation;
            context=(source=component,),
        ).rate
        rate = RateElement(
            process.formulation,
            (TracerOp{source}(), parameter_operand(rate_binding)),
        )
        loss = FluxSpec(process_id(named), source, rate, Weight{-1}())
        gain = FluxSpec(process_id(named), topology.destination_target, rate, Weight{1}())
        fluxes = (fluxes..., loss, gain)
    end
    return fluxes
end

function process_fluxes(
    named::NamedProcess{P},
    definition::NormalizedModelDefinition,
    layout::ComponentLayout,
    context::CommunityContext,
) where {P<:Remineralization}
    topology = realize_process_topology(named, layout, context)
    return process_fluxes(named, topology, definition)
end
