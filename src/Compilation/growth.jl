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

"""Resolved parameter names needed to compile one factorized growth process."""
struct GrowthParameterBinding{K,C,S,T}
    maximum_rate::Symbol
    half_saturation::K
    alpha::Symbol
    chlorophyll_to_carbon_ratio::C
    stoichiometry::S
    temperature::T
end

function _growth_stoichiometry_bindings(
    definition::NormalizedModelDefinition,
    named::NamedProcess,
    nutrients::Nutrients,
)
    isnothing(named.process.stoichiometry) && return NamedTuple()
    names = keys(nutrients.responses)
    values = Tuple(
        parameter_name(
            definition,
            _parameter_requirement(
                named,
                (:stoichiometry,),
                :ratio;
                qualifier=(currency=name,),
            ),
        ) for name in names
    )
    return NamedTuple{names}(values)
end

"""Resolve growth parameter names from normalized semantic requirement bindings."""
function GrowthParameterBinding(definition::NormalizedModelDefinition, id::Symbol)
    hasproperty(definition.processes, id) || throw(
        ArgumentError("normalized model has no process :$id"),
    )
    named = getproperty(definition.processes, id)
    process = named.process
    process isa Growth || throw(ArgumentError("process :$id is not a Growth process"))
    light_name, light = _growth_light_factor(process)
    nutrient_name, nutrients = _growth_nutrient_factor(process)
    temperature = _temperature_factor_binding(definition, named)

    light_path = (:factors, light_name)
    maximum_rate = parameter_name(
        definition, _parameter_requirement(named, light_path, :maximum_rate)
    )
    alpha = parameter_name(definition, _parameter_requirement(named, light_path, :alpha))

    if light isa Light{Smith} && nutrients isa NutrientResponse{Monod}
        half_saturation = parameter_name(
            definition, _parameter_requirement(named, (:factors, nutrient_name), :K;
                qualifier=(resource=nutrients.resource,))
        )
        return GrowthParameterBinding(
            maximum_rate, half_saturation, alpha, nothing, NamedTuple(), temperature
        )
    end

    if light isa Light{Geider} && nutrients isa Nutrients{Liebig}
        response_names = keys(nutrients.responses)
        half_saturations = NamedTuple{response_names}(Tuple(
            parameter_name(
                definition,
                _parameter_requirement(
                    named,
                    (:factors, nutrient_name, :responses, response_name),
                    :K;
                    qualifier=(resource=getproperty(nutrients.responses, response_name).resource,),
                ),
            ) for response_name in response_names
        ))
        chlorophyll = parameter_name(
            definition,
            _parameter_requirement(named, light_path, :chlorophyll_to_carbon_ratio),
        )
        stoichiometry = _growth_stoichiometry_bindings(definition, named, nutrients)
        return GrowthParameterBinding(
            maximum_rate, half_saturations, alpha, chlorophyll, stoichiometry, temperature
        )
    end

    throw(ArgumentError(
        "unsupported growth factor combination $(typeof(light)) × $(typeof(nutrients))",
    ))
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
    topology::GrowthTopology,
    binding::GrowthParameterBinding,
    population_index::Int,
)
    process = named.process
    _, light = _growth_light_factor(process)
    _, nutrients = _growth_nutrient_factor(process)
    formulations = (formulation(light), formulation(nutrients))
    rate_factors = _rate_factors(binding.temperature)

    if light isa Light{Smith} && nutrients isa NutrientResponse{Monod}
        operands = (
            ClassOp{population_index}(),
            TracerOp{topology.resource_target}(),
            TracerOp{topology.light_driver}(),
            VecParamOp{binding.maximum_rate,population_index}(),
            VecParamOp{binding.half_saturation,population_index}(),
            VecParamOp{binding.alpha,population_index}(),
        )
        return RateElement(formulations, operands; factors=rate_factors)
    end

    if light isa Light{Geider} && nutrients isa Nutrients{Liebig}
        resource_ops = Tuple(TracerOp{target}() for target in values(topology.resource_target))
        half_saturation_ops = Tuple(
            VecParamOp{name,population_index}() for name in values(binding.half_saturation)
        )
        operands = (
            ClassOp{population_index}(),
            TupleOp(resource_ops),
            TracerOp{topology.light_driver}(),
            VecParamOp{binding.maximum_rate,population_index}(),
            TupleOp(half_saturation_ops),
            VecParamOp{binding.alpha,population_index}(),
            VecParamOp{binding.chlorophyll_to_carbon_ratio,population_index}(),
        )
        return RateElement(formulations, operands; factors=rate_factors)
    end

    throw(ArgumentError(
        "unsupported growth factor combination $(typeof(light)) × $(typeof(nutrients))",
    ))
end

function _growth_resource_fluxes(
    named::NamedProcess,
    topology::GrowthTopology,
    binding::GrowthParameterBinding,
    rate::RateElement,
    nutrients::NutrientResponse{Monod},
)
    return (FluxSpec(process_id(named), topology.resource_target, rate, Weight{-1}()),)
end

function _growth_resource_fluxes(
    named::NamedProcess,
    topology::GrowthTopology,
    binding::GrowthParameterBinding,
    rate::RateElement,
    nutrients::Nutrients{Liebig},
)
    isnothing(topology.source_target) && throw(
        ArgumentError("Geider/Liebig growth requires a canonical source component"),
    )
    fluxes = (
        FluxSpec(process_id(named), topology.source_target, rate, Weight{-1}()),
    )
    for currency in keys(topology.resource_target)
        target = getproperty(topology.resource_target, currency)
        scale = getproperty(binding.stoichiometry, currency)
        fluxes = (
            fluxes...,
            FluxSpec(
                process_id(named),
                target,
                rate,
                Weight{-1}((ScalarParamOp{scale}(),)),
            ),
        )
    end
    return fluxes
end

"""Derive biomass-gain and resource-loss fluxes for factorized growth."""
function process_fluxes(
    named::NamedProcess{P}, topology::GrowthTopology, binding::GrowthParameterBinding
) where {P<:Growth}
    length(topology.population_tracers) == length(topology.population_indices) || throw(
        ArgumentError("growth topology tracer and index counts must match"),
    )
    _, nutrients = _growth_nutrient_factor(named.process)
    fluxes = ()

    for i in eachindex(topology.population_tracers)
        tracer = topology.population_tracers[i]
        rate = _growth_rate(named, topology, binding, topology.population_indices[i])
        biomass = FluxSpec(process_id(named), tracer, rate, Weight{1}())
        resources = _growth_resource_fluxes(named, topology, binding, rate, nutrients)
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
    binding = GrowthParameterBinding(definition, process_id(named))
    return process_fluxes(named, topology, binding)
end
