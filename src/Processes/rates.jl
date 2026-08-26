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

@inline factor_value(
    ::NormalizedDroop, internal, reference, minimum_quota, maximum_quota
) = normalized_droop_limitation(internal, reference, minimum_quota, maximum_quota)

@inline factor_value(::Liebig, limitations::Tuple) = liebig_minimum(limitations)
@inline factor_value(::FrankTNorm, limitations::Tuple, sharpness) =
    frank_tnorm(limitations; sharpness)

@inline factor_value(::Q10, temperature, q10, reference_temperature) =
    q10_temperature_factor(temperature, q10, reference_temperature)

"""Evaluate the unmodified population-growth scale before sibling factors."""
@inline process_rate(::MultiplicativeFactors, biomass, maximum_rate) =
    maximum_rate * biomass

"""Evaluate one quota-regulated external nutrient uptake rate."""
@inline function process_rate(
    ::QuotaRegulatedMonod,
    resource,
    internal,
    reference,
    maximum_rate,
    half_saturation,
    minimum_quota,
    maximum_quota,
    hill,
)
    return maximum_rate * reference * monod_limitation(resource, half_saturation) *
           quota_uptake_regulation(
               internal, reference, minimum_quota, maximum_quota, hill
           )
end

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
