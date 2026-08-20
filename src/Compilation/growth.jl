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

_growth_resource_target(factor::NutrientResponse, layout::ComponentLayout) =
    _scalar_component_target(layout, factor.resource)

function _growth_resource_target(factor::Nutrients, layout::ComponentLayout)
    names = keys(factor.responses)
    targets = Tuple(
        _scalar_component_target(layout, response.resource)
        for response in values(factor.responses)
    )
    return NamedTuple{names}(targets)
end

"""Realized factorized growth topology over concrete population classes."""
struct GrowthTopology{TR,IX,R,S,L}
    population_tracers::TR
    population_indices::IX
    resource_target::R
    source_target::S
    layout::L
end

"""Resolve a named factorized growth process onto the current concrete class layout."""
function realize_process_topology(
    named::NamedProcess{P}, layout::ComponentLayout, context::CommunityContext
) where {P<:Growth}
    process = named.process
    _, nutrient_factor = _growth_resource_factor(process)
    population_tracers, population_indices = _realize_population_classes(
        named, process.populations, layout, context
    )
    resources = _growth_resource_target(nutrient_factor, layout)
    source = isnothing(process.source) ? nothing : _scalar_component_target(layout, process.source)
    return GrowthTopology(population_tracers, population_indices, resources, source, layout)
end

function _growth_scale_binding(
    definition::NormalizedModelDefinition, named::NamedProcess
)
    matches = ()
    for (name, factor) in pairs(factors(named))
        any(slot -> slot.name === :maximum_rate, parameter_slots(formulation(factor))) ||
            continue
        slots = parameter_slot_bindings(
            definition, named, (:factors, name), formulation(factor)
        )
        matches = (matches..., slots.maximum_rate)
    end
    length(matches) == 1 || throw(ArgumentError(
        "growth process :$(process_id(named)) must declare exactly one factor-owned maximum_rate slot",
    ))
    return only(matches)
end

function _growth_rate(
    named::NamedProcess,
    definition::NormalizedModelDefinition,
    topology::GrowthTopology,
    population_index::Int,
    scale_binding::ParameterBinding,
)
    axis_indices = (population=population_index,)
    rate_factors = _factor_elements(definition, named, topology.layout, axis_indices)
    operands = (
        ClassOp{population_index}(),
        parameter_operand(scale_binding, population_index),
    )
    return RateElement(formulation(named.process), operands; factors=rate_factors)
end

function _growth_resource_fluxes(
    named::NamedProcess,
    definition::NormalizedModelDefinition,
    topology::GrowthTopology,
    rate::RateElement,
    nutrients::NutrientResponse,
)
    return (FluxSpec(process_id(named), topology.resource_target, rate, Weight{-1}()),)
end

function _growth_resource_fluxes(
    named::NamedProcess,
    definition::NormalizedModelDefinition,
    topology::GrowthTopology,
    rate::RateElement,
    nutrients::Nutrients,
)
    isnothing(topology.source_target) && throw(ArgumentError(
        "multi-resource growth requires a canonical source component",
    ))
    fluxes = (FluxSpec(process_id(named), topology.source_target, rate, Weight{-1}()),)
    stoichiometry = named.process.stoichiometry
    isnothing(stoichiometry) && throw(ArgumentError(
        "multi-resource growth requires stoichiometry for resource routing",
    ))
    for currency in keys(topology.resource_target)
        target = getproperty(topology.resource_target, currency)
        ratio = parameter_slot_bindings(
            definition,
            named,
            (:stoichiometry,),
            stoichiometry;
            context=(currency=currency,),
        ).ratio
        fluxes = (
            fluxes...,
            FluxSpec(
                process_id(named),
                target,
                rate,
                Weight{-1}((parameter_operand(ratio),)),
            ),
        )
    end
    return fluxes
end

"""Derive biomass-gain and resource-loss fluxes for factorized growth."""
function process_fluxes(
    named::NamedProcess{P},
    topology::GrowthTopology,
    definition::NormalizedModelDefinition,
) where {P<:Growth}
    _, nutrients = _growth_resource_factor(named.process)
    scale_binding = _growth_scale_binding(definition, named)
    fluxes = ()

    for i in eachindex(topology.population_tracers)
        tracer = topology.population_tracers[i]
        rate = _growth_rate(
            named, definition, topology, topology.population_indices[i], scale_binding
        )
        biomass = FluxSpec(process_id(named), tracer, rate, Weight{1}())
        resources = _growth_resource_fluxes(named, definition, topology, rate, nutrients)
        fluxes = (fluxes..., biomass, resources...)
    end
    return fluxes
end

function process_fluxes(
    named::NamedProcess{P},
    definition::NormalizedModelDefinition,
    layout::ComponentLayout,
    context::CommunityContext,
) where {P<:Growth}
    topology = realize_process_topology(named, layout, context)
    return process_fluxes(named, topology, definition)
end
