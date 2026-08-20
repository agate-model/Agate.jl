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

"""Positive biomass effect from one realized growth-rate element."""
struct GrowthBiomassContribution{R,K,C,S,L,N,T} <: AbstractProcessContribution
    process::Symbol
    target::Symbol
    population_index::Int
    resource::R
    light_driver::Symbol
    maximum_rate_parameter::Symbol
    half_saturation_parameter::K
    alpha_parameter::Symbol
    chlorophyll_to_carbon_ratio_parameter::C
    stoichiometry_parameters::S
    light::L
    nutrients::N
    temperature_factor::T
end

"""Resource loss coupled to one realized growth-rate element."""
struct GrowthResourceLossContribution{R,K,C,S,L,N,T,Q} <: AbstractProcessContribution
    process::Symbol
    target::Symbol
    population::Symbol
    population_index::Int
    resource::R
    light_driver::Symbol
    maximum_rate_parameter::Symbol
    half_saturation_parameter::K
    alpha_parameter::Symbol
    chlorophyll_to_carbon_ratio_parameter::C
    stoichiometry_parameters::S
    light::L
    nutrients::N
    temperature_factor::T
    scale_parameter::Q
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

function _growth_contribution_fields(
    named::NamedProcess,
    topology::GrowthTopology,
    binding::GrowthParameterBinding,
    population_index::Int,
)
    process = named.process
    _, light = _growth_light_factor(process)
    _, nutrients = _growth_nutrient_factor(process)
    return (
        process_id(named),
        population_index,
        topology.resource_target,
        topology.light_driver,
        binding.maximum_rate,
        binding.half_saturation,
        binding.alpha,
        binding.chlorophyll_to_carbon_ratio,
        binding.stoichiometry,
        light,
        nutrients,
        binding.temperature,
    )
end

function _growth_resource_contributions(
    named::NamedProcess,
    topology::GrowthTopology,
    binding::GrowthParameterBinding,
    tracer::Symbol,
    fields::Tuple,
    nutrients::NutrientResponse{Monod},
)
    resource = GrowthResourceLossContribution(
        fields[1], topology.resource_target, tracer, fields[2:end]..., nothing
    )
    return (resource,)
end

function _growth_resource_contributions(
    named::NamedProcess,
    topology::GrowthTopology,
    binding::GrowthParameterBinding,
    tracer::Symbol,
    fields::Tuple,
    nutrients::Nutrients{Liebig},
)
    isnothing(topology.source_target) && throw(
        ArgumentError("Geider/Liebig growth requires a canonical source component"),
    )
    source = GrowthResourceLossContribution(
        fields[1], topology.source_target, tracer, fields[2:end]..., nothing
    )
    resources = ()
    for currency in keys(topology.resource_target)
        target = getproperty(topology.resource_target, currency)
        scale = getproperty(binding.stoichiometry, currency)
        loss = GrowthResourceLossContribution(
            fields[1], target, tracer, fields[2:end]..., scale
        )
        resources = (resources..., loss)
    end
    return (source, resources...)
end

"""Derive coupled biomass-gain and resource-loss contributions for factorized growth."""
function process_contributions(
    named::NamedProcess{P}, topology::GrowthTopology, binding::GrowthParameterBinding
) where {P<:Growth}
    length(topology.population_tracers) == length(topology.population_indices) || throw(
        ArgumentError("growth topology tracer and index counts must match"),
    )
    _, nutrients = _growth_nutrient_factor(named.process)
    contributions = ()

    for i in eachindex(topology.population_tracers)
        tracer = topology.population_tracers[i]
        population_index = topology.population_indices[i]
        fields = _growth_contribution_fields(named, topology, binding, population_index)
        biomass = GrowthBiomassContribution(fields[1], tracer, fields[2:end]...)
        resources = _growth_resource_contributions(
            named, topology, binding, tracer, fields, nutrients
        )
        contributions = (contributions..., biomass, resources...)
    end
    return contributions
end

function process_contributions(
    named::NamedProcess{P},
    definition::NormalizedModelDefinition,
    layout::ComponentLayout,
    context::CommunityContext,
) where {P<:Growth}
    topology = realize_process_topology(named, layout, context)
    binding = GrowthParameterBinding(definition, process_id(named))
    return process_contributions(named, topology, binding)
end

@inline _growth_scale(bgc, ::Val{nothing}) = one(eltype(bgc.parameters.maximum_growth_rate))
@inline _growth_scale(bgc, ::Val{Q}) where {Q} = getproperty(bgc.parameters, Q)
@inline _growth_effect(::Val{:gain}, rate, scale) = rate
@inline _growth_effect(::Val{:loss}, rate, scale) = -scale * rate

struct SmithMonodGrowthTerm{I,R,D,M,K,A,Q,Effect,T}
    light_formulation::Smith
    nutrient_formulation::Monod
    temperature::T
end

@inline function (term::SmithMonodGrowthTerm{I,R,D,M,K,A,Q,Effect,T})(
    bgc, args
) where {I,R,D,M,K,A,Q,Effect,T}
    biomass = bgc.tracers.plankton(args, I)
    resource = getproperty(bgc.tracers, R)(args)
    light = getproperty(bgc.tracers, D)(args)
    maximum_rate = @inbounds getproperty(bgc.parameters, M)[I]
    half_saturation = @inbounds getproperty(bgc.parameters, K)[I]
    alpha = @inbounds getproperty(bgc.parameters, A)[I]
    rate = process_rate(
        term.light_formulation,
        term.nutrient_formulation,
        biomass,
        resource,
        light,
        maximum_rate,
        half_saturation,
        alpha,
    )
    rate = _apply_factor(term.temperature, bgc, args, rate)
    return _growth_effect(Val(Effect), rate, _growth_scale(bgc, Val(Q)))
end

@generated function _growth_tracer_values(bgc, args, ::Val{Names}) where {Names}
    expressions = [:(getproperty(bgc.tracers, $(QuoteNode(name)))(args)) for name in Names]
    return Expr(:tuple, expressions...)
end

@generated function _growth_parameter_values(bgc, ::Val{Names}, ::Val{I}) where {Names,I}
    expressions = [
        :(@inbounds getproperty(bgc.parameters, $(QuoteNode(name)))[$I]) for name in Names
    ]
    return Expr(:tuple, expressions...)
end

struct GeiderLiebigGrowthTerm{I,R,D,M,K,A,C,Q,Effect,T}
    light_formulation::Geider
    nutrient_formulation::Liebig
    temperature::T
end

@inline function (term::GeiderLiebigGrowthTerm{I,R,D,M,K,A,C,Q,Effect,T})(
    bgc, args
) where {I,R,D,M,K,A,C,Q,Effect,T}
    biomass = bgc.tracers.plankton(args, I)
    resources = _growth_tracer_values(bgc, args, Val(R))
    light = getproperty(bgc.tracers, D)(args)
    maximum_rate = @inbounds getproperty(bgc.parameters, M)[I]
    half_saturations = _growth_parameter_values(bgc, Val(K), Val(I))
    alpha = @inbounds getproperty(bgc.parameters, A)[I]
    chlorophyll = @inbounds getproperty(bgc.parameters, C)[I]
    rate = process_rate(
        term.light_formulation,
        term.nutrient_formulation,
        biomass,
        resources,
        light,
        maximum_rate,
        half_saturations,
        alpha,
        chlorophyll,
    )
    rate = _apply_factor(term.temperature, bgc, args, rate)
    return _growth_effect(Val(Effect), rate, _growth_scale(bgc, Val(Q)))
end

function _growth_term(contribution, ::Val{Effect}) where {Effect}
    light = contribution.light
    nutrients = contribution.nutrients
    scale = hasproperty(contribution, :scale_parameter) ? contribution.scale_parameter : nothing
    temperature = _lower_factor(contribution.temperature_factor)

    if light isa Light{Smith} && nutrients isa NutrientResponse{Monod}
        return SmithMonodGrowthTerm{
            contribution.population_index,
            contribution.resource,
            contribution.light_driver,
            contribution.maximum_rate_parameter,
            contribution.half_saturation_parameter,
            contribution.alpha_parameter,
            scale,
            Effect,
            typeof(temperature),
        }(formulation(light), formulation(nutrients), temperature)
    end

    if light isa Light{Geider} && nutrients isa Nutrients{Liebig}
        resource_names = Tuple(values(contribution.resource))
        half_saturation_names = Tuple(values(contribution.half_saturation_parameter))
        return GeiderLiebigGrowthTerm{
            contribution.population_index,
            resource_names,
            contribution.light_driver,
            contribution.maximum_rate_parameter,
            half_saturation_names,
            contribution.alpha_parameter,
            contribution.chlorophyll_to_carbon_ratio_parameter,
            scale,
            Effect,
            typeof(temperature),
        }(formulation(light), formulation(nutrients), temperature)
    end

    throw(ArgumentError(
        "unsupported growth contribution factor combination " *
        "$(typeof(light)) × $(typeof(nutrients))",
    ))
end

_lower_contribution(contribution::GrowthBiomassContribution) =
    _growth_term(contribution, Val(:gain))
_lower_contribution(contribution::GrowthResourceLossContribution) =
    _growth_term(contribution, Val(:loss))
