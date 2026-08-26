_growth_resource_factor(named::NamedProcess) =
    getproperty(named.process.factors, named.facts.routing.factor)

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
    name = named.facts.maximum_rate_factor
    factor = getproperty(factors(named), name)
    return parameter_slot_bindings(
        definition, named, (:factors, name), factor
    ).maximum_rate
end

function _growth_rate(
    named::NamedProcess,
    definition::NormalizedModelDefinition,
    layout::ModelLayout,
    plan::ParameterPlan,
    population_axis::Int,
    population_tracer::Symbol,
    scale_binding::ParameterBinding,
)
    axis_positions = (population=_axis_position(population_axis),)
    rate_factors = _factor_elements(definition, named, layout, plan, axis_positions)
    operands = (
        TracerOp{population_tracer}(),
        parameter_operand(scale_binding, plan, axis_positions),
    )
    return RateElement(formulation(named.process), operands; factors=rate_factors)
end

function _growth_resource_fluxes(
    named::NamedProcess,
    definition::NormalizedModelDefinition,
    plan::ParameterPlan,
    resource_target,
    source_target,
    rate::RateElement,
    ::NutrientResponse,
)
    return (FluxSpec(process_id(named), resource_target, rate, Weight{-1}()),)
end

function _growth_resource_fluxes(
    named::NamedProcess,
    definition::NormalizedModelDefinition,
    plan::ParameterPlan,
    resource_target,
    source_target,
    rate::RateElement,
    ::Nutrients,
)
    fluxes = Any[FluxSpec(process_id(named), source_target, rate, Weight{-1}())]
    stoichiometry = named.facts.routing.stoichiometry
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
                Weight{-1}((parameter_operand(ratio, plan),)),
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
    plan::ParameterPlan,
) where {P<:Growth}
    process = named.process
    nutrients = _growth_resource_factor(named)
    population_tracers = _realize_normalized_population_states(
        named.facts.population_states, layout
    )
    resource_target = _growth_resource_target(nutrients, layout)
    source_target = named.facts.routing.mode === :single_resource ? nothing :
                    _scalar_component_target(layout, named.facts.routing.source)
    scale_binding = _growth_scale_binding(definition, named)
    fluxes = Any[]

    for population_axis in eachindex(population_tracers)
        rate = _growth_rate(
            named, definition, layout, plan, population_axis,
            population_tracers[population_axis], scale_binding,
        )
        biomass = FluxSpec(
            process_id(named), population_tracers[population_axis], rate, Weight{1}()
        )
        resources = _growth_resource_fluxes(
            named, definition, plan, resource_target, source_target, rate, nutrients
        )
        push!(fluxes, biomass)
        append!(fluxes, resources)
    end
    return Tuple(fluxes)
end
