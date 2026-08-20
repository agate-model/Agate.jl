function _growth_factor(process::Growth, ::Type{F}, label::Symbol) where {F<:AbstractFactor}
    matches = ()
    for pair in pairs(process.factors)
        last(pair) isa F && (matches = (matches..., pair))
    end
    length(matches) == 1 || throw(
        ArgumentError("growth process must declare exactly one $label factor"),
    )
    pair = only(matches)
    return first(pair), last(pair)
end

"""Resolved parameter names needed to compile one Smith/Monod factorized growth process."""
struct GrowthParameterBinding
    maximum_rate::Symbol
    half_saturation::Symbol
    alpha::Symbol
end

"""Resolve growth parameter names from normalized semantic requirement bindings."""
function GrowthParameterBinding(definition::NormalizedModelDefinition, id::Symbol)
    hasproperty(definition.processes, id) || throw(
        ArgumentError("normalized model has no process :$id"),
    )
    named = getproperty(definition.processes, id)
    process = named.process
    process isa Growth || throw(ArgumentError("process :$id is not a Growth process"))
    light_name, _ = _growth_factor(process, Light{Smith}, :Smith_light)
    nutrient_name, _ = _growth_factor(process, NutrientResponse{Monod}, :Monod_nutrient)

    light_path = (:factors, light_name)
    nutrient_path = (:factors, nutrient_name)
    maximum_rate = parameter_name(
        definition, _parameter_requirement(named, light_path, :maximum_rate)
    )
    alpha = parameter_name(definition, _parameter_requirement(named, light_path, :alpha))
    half_saturation = parameter_name(
        definition, _parameter_requirement(named, nutrient_path, :K)
    )
    return GrowthParameterBinding(maximum_rate, half_saturation, alpha)
end

"""Realized Smith/Monod growth topology over concrete population classes."""
struct GrowthTopology{TR,IX,R,D}
    population_tracers::TR
    population_indices::IX
    resource_target::R
    light_driver::D
end

"""Positive biomass effect from one realized growth-rate element."""
struct GrowthBiomassContribution{L,N} <: AbstractProcessContribution
    process::Symbol
    target::Symbol
    population_index::Int
    resource::Symbol
    light_driver::Symbol
    maximum_rate_parameter::Symbol
    half_saturation_parameter::Symbol
    alpha_parameter::Symbol
    light::L
    nutrients::N
end

"""Resource loss coupled to one realized growth-rate element."""
struct GrowthResourceLossContribution{L,N} <: AbstractProcessContribution
    process::Symbol
    target::Symbol
    population::Symbol
    population_index::Int
    resource::Symbol
    light_driver::Symbol
    maximum_rate_parameter::Symbol
    half_saturation_parameter::Symbol
    alpha_parameter::Symbol
    light::L
    nutrients::N
end

"""Resolve a named Smith/Monod growth process onto the current concrete class layout."""
function realize_process_topology(
    named::NamedProcess{P}, layout::ComponentLayout, context::CommunityContext
) where {P<:Growth}
    process = named.process
    _, light = _growth_factor(process, Light{Smith}, :Smith_light)
    _, nutrients = _growth_factor(process, NutrientResponse{Monod}, :Monod_nutrient)

    population_tracers, population_indices = _realize_population_classes(
        named, process.populations, layout, context
    )
    resource = _scalar_component_target(layout, nutrients.resource)
    return GrowthTopology(
        population_tracers, population_indices, resource, light.driver
    )
end

"""Derive coupled biomass-gain and nutrient-loss contributions for Smith/Monod growth."""
function process_contributions(
    named::NamedProcess{P}, topology::GrowthTopology, binding::GrowthParameterBinding
) where {P<:Growth}
    length(topology.population_tracers) == length(topology.population_indices) || throw(
        ArgumentError("growth topology tracer and index counts must match"),
    )
    process = named.process
    _, light = _growth_factor(process, Light{Smith}, :Smith_light)
    _, nutrients = _growth_factor(process, NutrientResponse{Monod}, :Monod_nutrient)
    contributions = ()

    for i in eachindex(topology.population_tracers)
        tracer = topology.population_tracers[i]
        population_index = topology.population_indices[i]
        fields = (
            process_id(named),
            population_index,
            topology.resource_target,
            topology.light_driver,
            binding.maximum_rate,
            binding.half_saturation,
            binding.alpha,
            light,
            nutrients,
        )
        biomass = GrowthBiomassContribution(fields[1], tracer, fields[2:end]...)
        resource = GrowthResourceLossContribution(
            fields[1], topology.resource_target, tracer, fields[2:end]...
        )
        contributions = (contributions..., biomass, resource)
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

struct GrowthTerm{L,N,I,R,D,M,K,A,Effect}
    light_formulation::L
    nutrient_formulation::N
end

@inline _growth_effect(::Val{:gain}, rate) = rate
@inline _growth_effect(::Val{:loss}, rate) = -rate

@inline function (term::GrowthTerm{L,N,I,R,D,M,K,A,Effect})(
    bgc, args
) where {L,N,I,R,D,M,K,A,Effect}
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
    return _growth_effect(Val(Effect), rate)
end

function _growth_term(contribution, ::Val{Effect}) where {Effect}
    light_formulation = formulation(contribution.light)
    nutrient_formulation = formulation(contribution.nutrients)
    L = typeof(light_formulation)
    N = typeof(nutrient_formulation)
    return GrowthTerm{
        L,
        N,
        contribution.population_index,
        contribution.resource,
        contribution.light_driver,
        contribution.maximum_rate_parameter,
        contribution.half_saturation_parameter,
        contribution.alpha_parameter,
        Effect,
    }(light_formulation, nutrient_formulation)
end

_lower_contribution(contribution::GrowthBiomassContribution) =
    _growth_term(contribution, Val(:gain))
_lower_contribution(contribution::GrowthResourceLossContribution) =
    _growth_term(contribution, Val(:loss))
