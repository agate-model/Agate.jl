"""Resolved parameter name needed to compile one linear remineralization process."""
struct RemineralizationParameterBinding{R}
    rate::R
end

function RemineralizationParameterBinding(
    definition::NormalizedModelDefinition, id::Symbol
)
    hasproperty(definition.processes, id) || throw(
        ArgumentError("normalized model has no process :$id"),
    )
    named = getproperty(definition.processes, id)
    process = named.process
    process isa Remineralization || throw(
        ArgumentError("process :$id is not a Remineralization process"),
    )
    process.formulation isa LinearRemineralization || throw(
        ArgumentError("unsupported remineralization formulation $(typeof(process.formulation))"),
    )
    sources = process.sources
    rates = Tuple(
        parameter_name(
            definition,
            _parameter_requirement(named, (), :rate; qualifier=(source=source,)),
        ) for source in sources
    )
    rate = length(rates) == 1 ? only(rates) : NamedTuple{sources}(rates)
    return RemineralizationParameterBinding(rate)
end

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
    binding::RemineralizationParameterBinding,
) where {P<:Remineralization}
    process = named.process
    length(topology.source_tracers) == length(process.sources) || throw(
        ArgumentError("remineralization topology source count must match process sources"),
    )
    fluxes = ()
    for i in eachindex(topology.source_tracers)
        source = topology.source_tracers[i]
        component = process.sources[i]
        rate_parameter = binding.rate isa Symbol ? binding.rate : getproperty(binding.rate, component)
        rate = RateElement(
            process.formulation,
            (TracerOp{source}(), ScalarParamOp{rate_parameter}()),
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
    binding = RemineralizationParameterBinding(definition, process_id(named))
    return process_fluxes(named, topology, binding)
end
