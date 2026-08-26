function _mortality_slots(
    definition::NormalizedModelDefinition, named::NamedProcess, population::Symbol
)
    return parameter_slot_bindings(
        definition, named, (), named.process; context=(population=population,)
    )
end

function _mortality_rate(
    formulation,
    rate_binding::ParameterBinding,
    layout::ModelLayout,
    population_axis::Int,
    population_index::Int,
    population_tracer::Symbol,
)
    axis_positions = (population=_axis_position(population_axis, population_index),)
    operands = (
        TracerOp{population_tracer}(),
        parameter_operand(rate_binding, layout, axis_positions),
    )
    return RateElement(formulation, operands)
end

function process_fluxes(
    named::NamedProcess{P},
    definition::NormalizedModelDefinition,
    layout::ModelLayout,
) where {P<:Mortality}
    process = named.process
    fluxes = Any[]

    for reference in named.facts.population_states
        population = reference.population
        population_tracers, population_indices = _realize_normalized_population_states(
            (reference,), layout
        )
        slots = _mortality_slots(definition, named, population)
        for population_axis in eachindex(population_tracers)
            population_index = population_indices[population_axis]
            rate = _mortality_rate(
                formulation(process), slots.rate, layout, population_axis, population_index,
                population_tracers[population_axis],
            )
            push!(
                fluxes,
                FluxSpec(
                    process_id(named), population_tracers[population_axis], rate, Weight{-1}()
                ),
            )
            if process.products isa Products
                append!(fluxes, _product_fluxes(named, definition, process.products, layout, rate))
            end
        end
    end
    return Tuple(fluxes)
end
