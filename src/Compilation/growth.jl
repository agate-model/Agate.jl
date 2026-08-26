function _growth_resource_factor(process::Growth)
    matches = Tuple(
        pair for pair in pairs(process.factors)
        if last(pair) isa Union{NutrientResponse,Nutrients}
    )
    length(matches) == 1 || throw(ArgumentError(
        "growth process must declare exactly one nutrient factor for resource routing",
    ))
    return only(matches)
end

_growth_resource_target(factor::NutrientResponse, layout::ModelLayout) =
    _scalar_component_target(layout, factor.resource)

function _growth_resource_target(factor::Nutrients, layout::ModelLayout)
    names = keys(factor.responses)
    targets = Tuple(
        _scalar_component_target(layout, response.resource)
        for response in values(factor.responses)
    )
    return NamedTuple{names}(targets)
end

function _growth_scale_binding(
    definition::NormalizedModelDefinition, named::NamedProcess
)
    matches = ParameterBinding[]
    for (name, factor) in pairs(factors(named))
        any(slot -> slot.name === :maximum_rate, parameter_slots(formulation(factor))) ||
            continue
        slots = parameter_slot_bindings(
            definition, named, (:factors, name), factor
        )
        push!(matches, slots.maximum_rate)
    end
    length(matches) == 1 || throw(ArgumentError(
        "growth process :$(process_id(named)) must declare exactly one factor-owned maximum_rate slot",
    ))
    return only(matches)
end

function _growth_rate(
    named::NamedProcess,
    definition::NormalizedModelDefinition,
    layout::ModelLayout,
    population_axis::Int,
    population_index::Int,
    population_tracer::Symbol,
    scale_binding::ParameterBinding,
)
    axis_positions = (population=_axis_position(population_axis, population_index),)
    rate_factors = _factor_elements(definition, named, layout, axis_positions)
    operands = (
        TracerOp{population_tracer}(),
        parameter_operand(scale_binding, layout, axis_positions),
    )
    return RateElement(formulation(named.process), operands; factors=rate_factors)
end

function _growth_resource_fluxes(
    named::NamedProcess,
    definition::NormalizedModelDefinition,
    resource_target,
    source_target,
    rate::RateElement,
    nutrients::NutrientResponse,
)
    return (FluxSpec(process_id(named), resource_target, rate, Weight{-1}()),)
end

function _growth_resource_fluxes(
    named::NamedProcess,
    definition::NormalizedModelDefinition,
    resource_target,
    source_target,
    rate::RateElement,
    nutrients::Nutrients,
)
    isnothing(source_target) && throw(ArgumentError(
        "multi-resource growth requires a canonical source component",
    ))
    fluxes = Any[FluxSpec(process_id(named), source_target, rate, Weight{-1}())]
    stoichiometry = named.process.stoichiometry
    isnothing(stoichiometry) && throw(ArgumentError(
        "multi-resource growth requires stoichiometry for resource routing",
    ))
    for currency in keys(resource_target)
        ratio = parameter_slot_bindings(
            definition,
            named,
            (:stoichiometry,),
            stoichiometry;
            context=(currency=currency,),
        ).ratio
        push!(
            fluxes,
            FluxSpec(
                process_id(named),
                getproperty(resource_target, currency),
                rate,
                Weight{-1}((parameter_operand(ratio),)),
            ),
        )
    end
    return Tuple(fluxes)
end

"""Derive biomass-gain and resource-loss fluxes for factorized growth."""
function process_fluxes(
    named::NamedProcess{P},
    definition::NormalizedModelDefinition,
    layout::ModelLayout,
) where {P<:Growth}
    process = named.process
    _, nutrients = _growth_resource_factor(process)
    population_tracers, population_indices = _realize_population_classes(
        named, process.populations, layout
    )
    resource_target = _growth_resource_target(nutrients, layout)
    source_target = isnothing(process.source) ? nothing :
                    _scalar_component_target(layout, process.source)
    scale_binding = _growth_scale_binding(definition, named)
    fluxes = Any[]

    for population_axis in eachindex(population_tracers)
        population_index = population_indices[population_axis]
        rate = _growth_rate(
            named,
            definition,
            layout,
            population_axis,
            population_index,
            population_tracers[population_axis],
            scale_binding,
        )
        biomass = FluxSpec(
            process_id(named), population_tracers[population_axis], rate, Weight{1}()
        )
        resources = _growth_resource_fluxes(
            named, definition, resource_target, source_target, rate, nutrients
        )
        push!(fluxes, biomass)
        append!(fluxes, resources)
    end
    return Tuple(fluxes)
end
