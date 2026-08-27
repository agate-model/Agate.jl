function _growth_rate(
    named::NamedProcess,
    context::CompileContext,
    population_axis::Int,
    population_class::Symbol,
    population_tracer::Symbol,
    scale_ref::Int,
)
    axis_positions = (population=_axis_position(population_axis, population_class),)
    rate_factors = _factor_elements(context, named, axis_positions)
    operands = (
        input_operand(context.layout, population_tracer),
        parameter_operand(scale_ref, context, axis_positions),
    )
    return RateElement(formulation(named.process), operands; factors=rate_factors)
end

function _growth_resource_fluxes(
    named::NamedProcess,
    context::CompileContext,
    rate::RateElement,
)
    facts = named.facts
    layout = context.layout
    fluxes = Any[
        FluxSpec(
            _scalar_component_target(layout, facts.reference_source),
            rate,
            Weight{-1}(),
        ),
    ]
    for (currency, resource) in pairs(facts.additional_resources)
        ratio_ref = getproperty(named.binding_refs.stoichiometry, currency).ratio
        push!(
            fluxes,
            FluxSpec(
                _scalar_component_target(layout, resource),
                rate,
                Weight{-1}((parameter_operand(ratio_ref, context),)),
            ),
        )
    end
    return Tuple(fluxes)
end

"""Derive biomass-gain and resource-loss fluxes for factorized growth."""
function process_fluxes(
    named::NamedProcess{P}, context::CompileContext
) where {P<:Growth}
    layout = context.layout
    population_tracers = _realize_population_states(
        named.facts.population_states, layout
    )
    population_classes = _realize_population_classes(
        named.facts.population_states, layout
    )
    scale_ref = named.binding_refs.process.maximum_rate
    fluxes = Any[]

    for population_axis in eachindex(population_tracers)
        rate = _growth_rate(
            named, context, population_axis, population_classes[population_axis],
            population_tracers[population_axis], scale_ref,
        )
        push!(
            fluxes,
            FluxSpec(population_tracers[population_axis], rate, Weight{1}()),
        )
        append!(fluxes, _growth_resource_fluxes(named, context, rate))
    end
    return Tuple(fluxes)
end
