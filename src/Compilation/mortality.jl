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
    plan::ParameterPlan,
    population_axis::Int,
    population_tracer::Symbol,
)
    axis_positions = (population=_axis_position(population_axis),)
    operands = (
        TracerOp{population_tracer}(),
        parameter_operand(rate_binding, plan, axis_positions),
    )
    return RateElement(formulation, operands)
end

function process_fluxes(
    named::NamedProcess{P},
    definition::NormalizedModelDefinition,
    layout::ModelLayout,
    plan::ParameterPlan,
) where {P<:Mortality}
    process = named.process
    fluxes = Any[]

    for reference in named.facts.population_states
        population = reference.population
        population_tracers = _realize_normalized_population_states((reference,), layout)
        slots = _mortality_slots(definition, named, population)
        for population_axis in eachindex(population_tracers)
            rate = _mortality_rate(
                formulation(process), slots.rate, plan, population_axis,
                population_tracers[population_axis],
            )
            push!(
                fluxes,
                FluxSpec(
                    process_id(named), population_tracers[population_axis], rate, Weight{-1}()
                ),
            )
            if process.products isa Products
                append!(fluxes, _product_fluxes(named, definition, process.products, layout, plan, rate))
            end
        end
    end
    return Tuple(fluxes)
end
