function _mortality_slots(
    definition::NormalizedModelDefinition, named::NamedProcess
)
    populations = named.process.populations
    context = length(populations) == 1 ? (population=only(populations),) : NamedTuple()
    return parameter_slot_bindings(
        definition, named, (), formulation(named.process); context
    )
end

function _mortality_rate(formulation, rate_binding::ParameterBinding, population_index::Int)
    operands = (
        ClassOp{population_index}(),
        parameter_operand(rate_binding, population_index),
    )
    return RateElement(formulation, operands)
end

function process_fluxes(
    named::NamedProcess{P},
    definition::NormalizedModelDefinition,
    layout::ComponentLayout,
    context::CommunityContext,
) where {P<:Mortality}
    process = named.process
    population_tracers, population_indices = _realize_population_classes(
        named, process.populations, layout, context
    )
    slots = _mortality_slots(definition, named)
    fluxes = ()

    for i in eachindex(population_tracers)
        rate = _mortality_rate(formulation(process), slots.rate, population_indices[i])
        fluxes = (
            fluxes...,
            FluxSpec(process_id(named), population_tracers[i], rate, Weight{-1}()),
        )
        if !isnothing(process.routing)
            fluxes = (
                fluxes...,
                _routing_fluxes(named, definition, process.routing, layout, rate)...,
            )
        end
    end
    return fluxes
end
