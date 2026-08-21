"""Evaluate the canonical scientific rate for a mortality formulation."""
@inline process_rate(::LinearMortality, biomass, coefficient) = linear_loss(biomass, coefficient)
@inline process_rate(::QuadraticMortality, biomass, coefficient) = quadratic_loss(biomass, coefficient)
@inline factor_value(::Smith, light, maximum_rate, alpha) =
    smith_light_limitation(light, alpha, maximum_rate)

@inline factor_value(
    ::Geider, light, maximum_rate, alpha, chlorophyll_to_carbon_ratio
) = geider_light_response(light, alpha, maximum_rate, chlorophyll_to_carbon_ratio)

@inline factor_value(::Monod, resource, half_saturation) =
    monod_limitation(resource, half_saturation)

@inline factor_value(::Liebig, limitations::Tuple) = liebig_minimum(limitations)
@inline factor_value(formulation::Frank, limitations::Tuple) =
    FrankTNorm(formulation.sharpness)(limitations)

@inline factor_value(::Q10, temperature, q10, reference_temperature) =
    q10_temperature_factor(temperature, q10, reference_temperature)

"""Evaluate the unmodified population-growth scale before sibling factors."""
@inline process_rate(::MultiplicativeFactors, biomass, maximum_rate) =
    maximum_rate * biomass

"""Evaluate one heterotrophic consumer-by-resource uptake rate."""
@inline function process_rate(
    ::HeterotrophicConsumption, resource, consumer, maximum_rate, half_saturation
)
    return maximum_rate * factor_value(Monod(), resource, half_saturation) * consumer
end

"""Evaluate one idealized consumer-by-resource grazing rate."""
@inline function process_rate(
    ::IdealizedGrazing,
    resource,
    consumer,
    maximum_rate,
    half_saturation,
    palatability,
)
    return palatability * idealized_predation_loss(
        resource, consumer, maximum_rate, half_saturation
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

"""Evaluate one linear source remineralization rate."""
@inline process_rate(::LinearRemineralization, source, coefficient) =
    linear_remineralization(source, coefficient)
