"""Resolved parameter names for one heterotrophic consumption process."""
struct ConsumptionParameterBinding{T}
    maximum_rate::Symbol
    half_saturation::Symbol
    assimilation::Symbol
    temperature::T
end

function ConsumptionParameterBinding(definition::NormalizedModelDefinition, id::Symbol)
    hasproperty(definition.processes, id) || throw(
        ArgumentError("normalized model has no process :$id"),
    )
    named = getproperty(definition.processes, id)
    process = named.process
    process isa Consumption || throw(ArgumentError("process :$id is not a Consumption process"))
    process.formulation isa HeterotrophicConsumption || throw(
        ArgumentError("unsupported consumption formulation $(typeof(process.formulation))"),
    )
    return ConsumptionParameterBinding(
        parameter_name(definition, _parameter_requirement(named, (), :maximum_rate)),
        parameter_name(definition, _parameter_requirement(named, (), :half_saturation)),
        parameter_name(definition, _parameter_requirement(named, (), :assimilation)),
        _temperature_factor_binding(definition, named),
    )
end

"""Realized consumer-by-resource topology for heterotrophic consumption."""
struct ConsumptionTopology{CT,CI,RT,D}
    consumer_tracers::CT
    consumer_indices::CI
    resource_tracers::RT
    unassimilated_target::D
end

function _consumption_resource_tracers(
    named::NamedProcess, resources::Tuple, layout::ComponentLayout
)
    tracers = ()
    for resource in resources
        hasproperty(layout.component_tracers, resource) || throw(
            ArgumentError(
                "process :$(process_id(named)) references unrealized resource :$resource"
            ),
        )
        tracers = (tracers..., getproperty(layout.component_tracers, resource)...)
    end
    return tracers
end

function realize_process_topology(
    named::NamedProcess{P}, layout::ComponentLayout, context::CommunityContext
) where {P<:Consumption}
    process = named.process
    process.formulation isa HeterotrophicConsumption || throw(
        ArgumentError("unsupported consumption formulation $(typeof(process.formulation))"),
    )
    consumer_tracers, consumer_indices = _realize_population_classes(
        named, process.consumers, layout, context
    )
    resources = _consumption_resource_tracers(named, process.resources, layout)
    destination = isnothing(process.unassimilated_destination) ? nothing :
                  _scalar_component_target(layout, process.unassimilated_destination)
    return ConsumptionTopology(consumer_tracers, consumer_indices, resources, destination)
end

function _consumption_rate(
    formulation,
    binding::ConsumptionParameterBinding,
    consumer_index::Int,
    consumer_axis::Int,
    resource::Symbol,
    resource_axis::Int,
)
    operands = (
        TracerOp{resource}(),
        ClassOp{consumer_index}(),
        VecParamOp{binding.maximum_rate,consumer_axis}(),
        VecParamOp{binding.half_saturation,resource_axis}(),
    )
    return RateElement(formulation, operands; factors=_rate_factors(binding.temperature))
end

function process_fluxes(
    named::NamedProcess{P},
    topology::ConsumptionTopology,
    binding::ConsumptionParameterBinding,
) where {P<:Consumption}
    fluxes = ()
    for consumer_axis in eachindex(topology.consumer_tracers)
        consumer = topology.consumer_tracers[consumer_axis]
        consumer_index = topology.consumer_indices[consumer_axis]
        for resource_axis in eachindex(topology.resource_tracers)
            resource = topology.resource_tracers[resource_axis]
            rate = _consumption_rate(
                formulation(named.process),
                binding,
                consumer_index,
                consumer_axis,
                resource,
                resource_axis,
            )
            assimilation = MatParamOp{binding.assimilation,consumer_axis,resource_axis}()
            loss = FluxSpec(process_id(named), resource, rate, Weight{-1}())
            gain = FluxSpec(process_id(named), consumer, rate, Weight{1}((assimilation,)))
            fluxes = (fluxes..., loss, gain)
            if !isnothing(topology.unassimilated_target)
                unassimilated = FluxSpec(
                    process_id(named),
                    topology.unassimilated_target,
                    rate,
                    Weight{1}((ComplementOp(assimilation),)),
                )
                fluxes = (fluxes..., unassimilated)
            end
        end
    end
    return fluxes
end

function process_fluxes(
    named::NamedProcess{P},
    definition::NormalizedModelDefinition,
    layout::ComponentLayout,
    context::CommunityContext,
) where {P<:Consumption}
    topology = realize_process_topology(named, layout, context)
    binding = ConsumptionParameterBinding(definition, process_id(named))
    return process_fluxes(named, topology, binding)
end
