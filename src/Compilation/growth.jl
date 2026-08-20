function _growth_light_factor(process::Growth)
    matches = ()
    for pair in pairs(process.factors)
        last(pair) isa Light && (matches = (matches..., pair))
    end
    length(matches) == 1 || throw(
        ArgumentError("growth process must declare exactly one light factor"),
    )
    pair = only(matches)
    return first(pair), last(pair)
end

function _growth_nutrient_factor(process::Growth)
    matches = ()
    for pair in pairs(process.factors)
        last(pair) isa Union{NutrientResponse,Nutrients} && (matches = (matches..., pair))
    end
    length(matches) == 1 || throw(
        ArgumentError("growth process must declare exactly one nutrient factor"),
    )
    pair = only(matches)
    return first(pair), last(pair)
end

"""Realized factorized growth topology over concrete population classes."""
struct GrowthTopology{TR,IX,R,S,D}
    population_tracers::TR
    population_indices::IX
    resource_target::R
    source_target::S
    light_driver::D
end

function _growth_resources(nutrients::NutrientResponse{Monod}, layout::ComponentLayout)
    return _scalar_component_target(layout, nutrients.resource)
end

function _growth_resources(nutrients::Nutrients{Liebig}, layout::ComponentLayout)
    names = keys(nutrients.responses)
    values = Tuple(
        _scalar_component_target(layout, getproperty(nutrients.responses, name).resource)
        for name in names
    )
    return NamedTuple{names}(values)
end

"""Resolve a named factorized growth process onto the current concrete class layout."""
function realize_process_topology(
    named::NamedProcess{P}, layout::ComponentLayout, context::CommunityContext
) where {P<:Growth}
    process = named.process
    _, light = _growth_light_factor(process)
    _, nutrients = _growth_nutrient_factor(process)

    population_tracers, population_indices = _realize_population_classes(
        named, process.populations, layout, context
    )
    resources = _growth_resources(nutrients, layout)
    source = isnothing(process.source) ? nothing : _scalar_component_target(layout, process.source)
    return GrowthTopology(
        population_tracers, population_indices, resources, source, light.driver
    )
end

function _growth_rate(
    named::NamedProcess,
    definition::NormalizedModelDefinition,
    topology::GrowthTopology,
    population_index::Int,
)
    process = named.process
    light_name, light = _growth_light_factor(process)
    nutrient_name, nutrients = _growth_nutrient_factor(process)
    formulations = (formulation(light), formulation(nutrients))
    light_slots = parameter_slot_bindings(
        definition, named, (:factors, light_name), light.formulation
    )
    rate_factors = _rate_factors(definition, named)

    if light isa Light{Smith} && nutrients isa NutrientResponse{Monod}
        nutrient_slots = parameter_slot_bindings(
            definition,
            named,
            (:factors, nutrient_name),
            nutrients.formulation;
            context=(resource=nutrients.resource,),
        )
        operands = (
            ClassOp{population_index}(),
            TracerOp{topology.resource_target}(),
            TracerOp{topology.light_driver}(),
            parameter_operand(light_slots.maximum_rate, population_index),
            parameter_operand(nutrient_slots.K, population_index),
            parameter_operand(light_slots.alpha, population_index),
        )
        return RateElement(formulations, operands; factors=rate_factors)
    end

    if light isa Light{Geider} && nutrients isa Nutrients{Liebig}
        resource_ops = Tuple(TracerOp{target}() for target in values(topology.resource_target))
        half_saturation_ops = Tuple(
            parameter_operand(
                parameter_slot_bindings(
                    definition,
                    named,
                    (:factors, nutrient_name, :responses, response_name),
                    response.formulation;
                    context=(resource=response.resource,),
                ).K,
                population_index,
            )
            for (response_name, response) in pairs(nutrients.responses)
        )
        operands = (
            ClassOp{population_index}(),
            TupleOp(resource_ops),
            TracerOp{topology.light_driver}(),
            parameter_operand(light_slots.maximum_rate, population_index),
            TupleOp(half_saturation_ops),
            parameter_operand(light_slots.alpha, population_index),
            parameter_operand(light_slots.chlorophyll_to_carbon_ratio, population_index),
        )
        return RateElement(formulations, operands; factors=rate_factors)
    end

    throw(ArgumentError(
        "unsupported growth factor combination $(typeof(light)) × $(typeof(nutrients))",
    ))
end

function _growth_resource_fluxes(
    named::NamedProcess,
    definition::NormalizedModelDefinition,
    topology::GrowthTopology,
    rate::RateElement,
    nutrients::NutrientResponse{Monod},
)
    return (FluxSpec(process_id(named), topology.resource_target, rate, Weight{-1}()),)
end

function _growth_resource_fluxes(
    named::NamedProcess,
    definition::NormalizedModelDefinition,
    topology::GrowthTopology,
    rate::RateElement,
    nutrients::Nutrients{Liebig},
)
    isnothing(topology.source_target) && throw(
        ArgumentError("Geider/Liebig growth requires a canonical source component"),
    )
    fluxes = (
        FluxSpec(process_id(named), topology.source_target, rate, Weight{-1}()),
    )
    stoichiometry = named.process.stoichiometry
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
    length(topology.population_tracers) == length(topology.population_indices) || throw(
        ArgumentError("growth topology tracer and index counts must match"),
    )
    _, nutrients = _growth_nutrient_factor(named.process)
    fluxes = ()

    for i in eachindex(topology.population_tracers)
        tracer = topology.population_tracers[i]
        rate = _growth_rate(named, definition, topology, topology.population_indices[i])
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
