"""Evaluate the canonical scientific rate for a mortality formulation."""
@inline process_rate(::LinearMortality, biomass, coefficient) = linear_loss(biomass, coefficient)
@inline process_rate(::QuadraticMortality, biomass, coefficient) = quadratic_loss(biomass, coefficient)
@inline process_rate(process::Mortality, biomass, coefficient) =
    process_rate(formulation(process), biomass, coefficient)

"""Evaluate Smith growth with one Monod nutrient response."""
@inline function process_rate(
    ::Smith,
    ::Monod,
    biomass,
    resource,
    light,
    maximum_rate,
    half_saturation,
    alpha,
)
    return smith_growth(
        (resource,), biomass, light, maximum_rate, (half_saturation,), alpha
    )
end

@inline function process_rate(
    process::Growth{Smith,L},
    biomass,
    resource,
    light,
    maximum_rate,
    half_saturation,
    alpha,
) where {L<:NutrientResponse{Monod}}
    return process_rate(
        process.formulation,
        process.limitation.formulation,
        biomass,
        resource,
        light,
        maximum_rate,
        half_saturation,
        alpha,
    )
end

"""Evaluate one preferential consumer-by-resource grazing rate."""
@inline function process_rate(
    ::PreferentialGrazing,
    resource,
    consumer,
    maximum_rate,
    half_saturation,
    palatability,
)
    return preferential_predation_loss(
        resource, consumer, maximum_rate, half_saturation, palatability
    )
end

@inline function process_rate(
    process::Grazing{PreferentialGrazing},
    resource,
    consumer,
    maximum_rate,
    half_saturation,
    palatability,
)
    return process_rate(
        process.formulation,
        resource,
        consumer,
        maximum_rate,
        half_saturation,
        palatability,
    )
end
