"""Resolved parameter names needed to compile one Smith/Monod growth process."""
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
    process.formulation isa Smith || throw(
        ArgumentError("unsupported growth formulation $(typeof(process.formulation))"),
    )
    process.limitation isa NutrientResponse{Monod} || throw(
        ArgumentError("growth compiler currently requires a Monod NutrientResponse"),
    )

    maximum_rate = parameter_name(
        definition, _parameter_requirement(named, (), :maximum_rate)
    )
    alpha = parameter_name(definition, _parameter_requirement(named, (), :alpha))
    half_saturation = parameter_name(
        definition, _parameter_requirement(named, (:limitation,), :K)
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
struct GrowthBiomassContribution{F,L} <: AbstractProcessContribution
    process::Symbol
    target::Symbol
    population_index::Int
    resource::Symbol
    light_driver::Symbol
    maximum_rate_parameter::Symbol
    half_saturation_parameter::Symbol
    alpha_parameter::Symbol
    formulation::F
    limitation::L
end

"""Resource loss coupled to one realized growth-rate element."""
struct GrowthResourceLossContribution{F,L} <: AbstractProcessContribution
    process::Symbol
    target::Symbol
    population::Symbol
    population_index::Int
    resource::Symbol
    light_driver::Symbol
    maximum_rate_parameter::Symbol
    half_saturation_parameter::Symbol
    alpha_parameter::Symbol
    formulation::F
    limitation::L
end

"""Resolve a named Smith/Monod growth process onto the current concrete class layout."""
function realize_process_topology(
    named::NamedProcess{P}, layout::ComponentLayout, context::CommunityContext
) where {P<:Growth}
    process = named.process
    process.formulation isa Smith || throw(
        ArgumentError("unsupported growth formulation $(typeof(process.formulation))"),
    )
    limitation = process.limitation
    limitation isa NutrientResponse{Monod} || throw(
        ArgumentError("growth compiler currently requires a Monod NutrientResponse"),
    )

    population_tracers, population_indices = _realize_population_classes(
        named, process.populations, layout, context
    )
    resource = _scalar_component_target(layout, limitation.resource)
    return GrowthTopology(
        population_tracers, population_indices, resource, process.light
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
    limitation = process.limitation
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
            process.formulation,
            limitation.formulation,
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

struct GrowthTerm{F,L,I,R,D,M,K,A,Effect}
    formulation::F
    limitation::L
end

@inline _growth_effect(::Val{:gain}, rate) = rate
@inline _growth_effect(::Val{:loss}, rate) = -rate

@inline function (term::GrowthTerm{F,L,I,R,D,M,K,A,Effect})(
    bgc, args
) where {F,L,I,R,D,M,K,A,Effect}
    biomass = bgc.tracers.plankton(args, I)
    resource = getproperty(bgc.tracers, R)(args)
    light = getproperty(bgc.tracers, D)(args)
    maximum_rate = @inbounds getproperty(bgc.parameters, M)[I]
    half_saturation = @inbounds getproperty(bgc.parameters, K)[I]
    alpha = @inbounds getproperty(bgc.parameters, A)[I]
    rate = process_rate(
        term.formulation,
        term.limitation,
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
    F = typeof(contribution.formulation)
    L = typeof(contribution.limitation)
    return GrowthTerm{
        F,
        L,
        contribution.population_index,
        contribution.resource,
        contribution.light_driver,
        contribution.maximum_rate_parameter,
        contribution.half_saturation_parameter,
        contribution.alpha_parameter,
        Effect,
    }(contribution.formulation, contribution.limitation)
end

_lower_contribution(contribution::GrowthBiomassContribution) =
    _growth_term(contribution, Val(:gain))
_lower_contribution(contribution::GrowthResourceLossContribution) =
    _growth_term(contribution, Val(:loss))
