"""Evaluate the canonical scientific rate for a mortality formulation."""
@inline process_rate(::LinearMortality, biomass, coefficient) = linear_loss(biomass, coefficient)
@inline process_rate(::QuadraticMortality, biomass, coefficient) = quadratic_loss(biomass, coefficient)
@inline process_rate(process::Mortality, biomass, coefficient) =
    process_rate(formulation(process), biomass, coefficient)

@inline factor_value(::Smith, light, maximum_rate, alpha) =
    smith_light_limitation(light, alpha, maximum_rate)

@inline factor_value(::Monod, resource, half_saturation) =
    monod_limitation(resource, half_saturation)

@inline factor_value(::Q10, temperature, q10, reference_temperature) =
    q10 ^ ((temperature - reference_temperature) / 10)

"""Evaluate growth by multiplying its Smith light and Monod nutrient factors."""
@inline function process_rate(
    light_formulation::Smith,
    nutrient_formulation::Monod,
    biomass,
    resource,
    light,
    maximum_rate,
    half_saturation,
    alpha,
)
    nutrient = factor_value(nutrient_formulation, resource, half_saturation)
    irradiance = factor_value(light_formulation, light, maximum_rate, alpha)
    return maximum_rate * nutrient * irradiance * biomass
end

"""Evaluate Geider growth with Liebig composition of Monod nutrient responses."""
@inline function process_rate(
    ::Geider,
    ::Liebig,
    biomass,
    resources::Tuple,
    light,
    maximum_rate,
    half_saturations::Tuple,
    alpha,
    chlorophyll_to_carbon_ratio,
)
    return geider_growth(
        resources,
        biomass,
        light,
        maximum_rate,
        half_saturations,
        alpha,
        chlorophyll_to_carbon_ratio,
    )
end

"""Evaluate one heterotrophic consumer-by-resource uptake rate."""
@inline function process_rate(
    ::HeterotrophicConsumption, resource, consumer, maximum_rate, half_saturation
)
    return maximum_rate * factor_value(Monod(), resource, half_saturation) * consumer
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

"""Evaluate one linear source remineralization rate."""
@inline process_rate(::LinearRemineralization, source, coefficient) =
    linear_remineralization(source, coefficient)
@inline process_rate(process::Remineralization{LinearRemineralization}, source, coefficient) =
    process_rate(process.formulation, source, coefficient)
