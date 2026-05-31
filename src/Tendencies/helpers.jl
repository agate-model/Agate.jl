@inline tracer_value(values, ::Val{name}) where {name} = getproperty(values, name)
@inline parameter_value(parameters, ::Val{name}) where {name} = getproperty(parameters, name)
@inline unity_from(x) = one(x)

@inline function stoichiometry_multiplier(parameters, ::Val{:one})
    return unity_from(parameters.maximum_growth_rate[1])
end
@inline stoichiometry_multiplier(parameters, ::Val{name}) where {name} = parameter_value(
    parameters, Val(name)
)

@inline function organic_fraction(parameters, ::Val{:DOM})
    f = parameters.DOM_POM_fractionation
    return one(f) - f
end
@inline organic_fraction(parameters, ::Val{:POM}) = parameters.DOM_POM_fractionation

@inline function half_saturation_parameter(parameters, ::NutrientCoupling{Tracer,HalfSaturation,Stoichiometry,Remineralization}) where {Tracer,HalfSaturation,Stoichiometry,Remineralization}
    return parameter_value(parameters, Val(HalfSaturation))
end

@inline function nutrient_resource(values, ::NutrientCoupling{Tracer,HalfSaturation,Stoichiometry,Remineralization}) where {Tracer,HalfSaturation,Stoichiometry,Remineralization}
    return tracer_value(values, Val(Tracer))
end

@inline nutrient_resources(values, nutrients::Tuple) = map(nutrient -> nutrient_resource(values, nutrient), nutrients)
@inline half_saturation_parameters(parameters, nutrients::Tuple) = map(nutrient -> half_saturation_parameter(parameters, nutrient), nutrients)

@inline function remineralization_sum(values, parameters, sources::Tuple)
    acc = zero(parameters.maximum_growth_rate[1])
    @inbounds for source in sources
        tracer_name, parameter_name = source
        tracer = tracer_value(values, Val(tracer_name))
        remin = parameter_value(parameters, Val(parameter_name))
        acc += linear_remineralization(tracer, remin)
    end
    return acc
end

@inline function growth_parameters(::Val{:smith_detritus}, parameters)
    return (parameters.alpha,)
end

@inline function growth_parameters(::Val{:geider_dom_pom}, parameters)
    return (parameters.photosynthetic_slope, parameters.chlorophyll_to_carbon_ratio)
end

@inline function plankton_growth(
    ::Val{:smith_detritus},
    ::Val{:liebig},
    resources::Tuple,
    P,
    PAR,
    parameters,
    half_saturations::Tuple,
    plankton_idx::Int,
)
    μmax = @inbounds parameters.maximum_growth_rate[plankton_idx]
    α = @inbounds parameters.alpha[plankton_idx]
    return smith_growth(resources, P, PAR, μmax, half_saturations, α)
end

@inline function plankton_growth(
    ::Val{:geider_dom_pom},
    ::Val{:liebig},
    resources::Tuple,
    P,
    PAR,
    parameters,
    half_saturations::Tuple,
    plankton_idx::Int,
)
    μmax = @inbounds parameters.maximum_growth_rate[plankton_idx]
    α = @inbounds parameters.photosynthetic_slope[plankton_idx]
    θc = @inbounds parameters.chlorophyll_to_carbon_ratio[plankton_idx]
    return geider_growth(resources, P, PAR, μmax, half_saturations, α, θc)
end
