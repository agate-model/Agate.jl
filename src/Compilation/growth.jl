_growth_resource_factor(named::NamedProcess) =
    getproperty(named.process.factors, named.facts.routing.factor)

_growth_resource_target(
    ::QuotaRouting, ::Nutrients, ::ModelLayout
) = nothing

_growth_resource_target(
    ::SingleResourceRouting, factor::NutrientResponse, layout::ModelLayout
) = _scalar_component_target(layout, factor.resource)

function _growth_resource_target(
    ::MultiResourceRouting, factor::Nutrients, layout::ModelLayout
)
    names = keys(factor.responses)
    targets = Tuple(
        _scalar_component_target(layout, response.resource)
        for response in values(factor.responses)
    )
    return NamedTuple{names}(targets)
end

_growth_source_target(::SingleResourceRouting, ::Growth, ::ModelLayout) = nothing
_growth_source_target(
    ::Union{QuotaRouting,MultiResourceRouting}, process::Growth, layout::ModelLayout
) = _scalar_component_target(layout, process.source)

function _growth_scale_binding(context::CompileContext, named::NamedProcess)
    name = named.facts.maximum_rate_factor
    factor = getproperty(factors(named), name)
    return parameter_slot_bindings(
        context.definition, named, (:factors, name), factor
    ).maximum_rate
end

function _growth_rate(
    named::NamedProcess,
    context::CompileContext,
    population_axis::Int,
    population_tracer::Symbol,
    scale_binding::ParameterBinding,
)
    axis_positions = (population=_axis_position(population_axis),)
    rate_factors = _factor_elements(context, named, axis_positions)
    operands = (
        input_operand(context.layout, population_tracer),
        parameter_operand(scale_binding, context.plan, axis_positions),
    )
    return RateElement(formulation(named.process), operands; factors=rate_factors)
end

function _growth_resource_fluxes(
    ::SingleResourceRouting,
    ::NamedProcess,
    ::CompileContext,
    resource_target,
    source_target,
    rate::RateElement,
    ::NutrientResponse,
)
    return (FluxSpec(resource_target, rate, Weight{-1}()),)
end

function _growth_resource_fluxes(
    ::QuotaRouting,
    ::NamedProcess,
    ::CompileContext,
    resource_target,
    source_target,
    rate::RateElement,
    ::Nutrients,
)
    return (FluxSpec(source_target, rate, Weight{-1}()),)
end

function _growth_resource_fluxes(
    route::MultiResourceRouting,
    named::NamedProcess,
    context::CompileContext,
    resource_target,
    source_target,
    rate::RateElement,
    ::Nutrients,
)
    fluxes = Any[FluxSpec(source_target, rate, Weight{-1}())]
    for currency in keys(resource_target)
        ratio = parameter_slot_bindings(
            context.definition,
            named,
            (:stoichiometry,),
            named.process.stoichiometry;
            context=(currency=currency,),
        ).ratio
        push!(
            fluxes,
            FluxSpec(
                getproperty(resource_target, currency),
                rate,
                Weight{-1}((parameter_operand(ratio, context.plan),)),
            ),
        )
    end
    return Tuple(fluxes)
end

"""Derive biomass-gain and resource-loss fluxes for factorized growth."""
function process_fluxes(
    named::NamedProcess{P}, context::CompileContext
) where {P<:Growth}
    route = named.facts.routing
    nutrients = _growth_resource_factor(named)
    layout = context.layout
    population_tracers = _realize_normalized_population_states(
        named.facts.population_states, layout
    )
    resource_target = _growth_resource_target(route, nutrients, layout)
    source_target = _growth_source_target(route, named.process, layout)
    scale_binding = _growth_scale_binding(context, named)
    fluxes = Any[]

    for population_axis in eachindex(population_tracers)
        rate = _growth_rate(
            named, context, population_axis,
            population_tracers[population_axis], scale_binding,
        )
        biomass = FluxSpec(
            population_tracers[population_axis], rate, Weight{1}()
        )
        resources = _growth_resource_fluxes(
            route, named, context, resource_target, source_target, rate, nutrients
        )
        push!(fluxes, biomass)
        append!(fluxes, resources)
    end
    return Tuple(fluxes)
end
