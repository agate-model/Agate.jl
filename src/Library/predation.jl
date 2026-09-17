"""Predation and grazing kernels."""
module Predation

export holling_type_ii, preferential_predation_loss

"""
    holling_type_ii(P, K)

Return the Holling (1959) type-II functional response ``P / (K + P)``.
The indeterminate `P == K == 0` case returns zero.
"""
@inline function holling_type_ii(P, K)
    K == zero(K) && P == zero(P) && return zero(P)
    return P / (K + P)
end

@inline _palatable_biomass(
    reference_inventories::Tuple{R}, palatabilities::Tuple{P}
) where {R,P} = first(reference_inventories) * first(palatabilities)

@inline function _palatable_biomass(
    reference_inventories::Tuple{R1,R2,Vararg{Any,N}},
    palatabilities::Tuple{P1,P2,Vararg{Any,N}},
) where {R1,R2,P1,P2,N}
    return first(reference_inventories) * first(palatabilities) + _palatable_biomass(
        Base.tail(reference_inventories), Base.tail(palatabilities)
    )
end

@inline _switching_weight_sum(
    reference_inventories::Tuple{R}, palatabilities::Tuple{P}, switching_exponent
) where {R,P} =
    (first(reference_inventories) * first(palatabilities))^switching_exponent

@inline function _switching_weight_sum(
    reference_inventories::Tuple{R1,R2,Vararg{Any,N}},
    palatabilities::Tuple{P1,P2,Vararg{Any,N}},
    switching_exponent,
) where {R1,R2,P1,P2,N}
    current = (first(reference_inventories) * first(palatabilities))^switching_exponent
    return current + _switching_weight_sum(
        Base.tail(reference_inventories), Base.tail(palatabilities), switching_exponent
    )
end

"""
    preferential_predation_loss(
        inventory, reference_inventory, consumer, maximum_grazing_rate,
        half_saturation, palatability, reference_inventories, palatabilities,
        switching_exponent
    )

Return the loss from one prey state when a consumer shares one maximum ingestion capacity across
all prey. Consumer-level saturation depends on total palatable reference biomass. A switching
exponent of one uses the algebraically simplified proportional-allocation path.
"""
@inline function preferential_predation_loss(
    inventory,
    reference_inventory,
    consumer,
    maximum_grazing_rate,
    half_saturation,
    palatability,
    reference_inventories::Tuple,
    palatabilities::Tuple,
    switching_exponent,
)
    palatable_biomass = _palatable_biomass(reference_inventories, palatabilities)
    if switching_exponent == one(switching_exponent)
        half_saturation == zero(half_saturation) && palatable_biomass == zero(palatable_biomass) &&
            return zero(maximum_grazing_rate * inventory * consumer)
        return maximum_grazing_rate * palatability * inventory /
               (half_saturation + palatable_biomass) * consumer
    end

    reference_inventory == zero(reference_inventory) &&
        return zero(maximum_grazing_rate * inventory * consumer)
    saturation = holling_type_ii(palatable_biomass, half_saturation)
    weights = _switching_weight_sum(
        reference_inventories, palatabilities, switching_exponent
    )
    weights == zero(weights) && return zero(maximum_grazing_rate * inventory * consumer)
    allocation = (palatability * reference_inventory)^switching_exponent / weights
    return maximum_grazing_rate * consumer * saturation * allocation *
           inventory / reference_inventory
end

@inline preferential_predation_loss(
    inventory, reference_inventory, consumer, maximum_grazing_rate, half_saturation, palatability
) = preferential_predation_loss(
    inventory,
    reference_inventory,
    consumer,
    maximum_grazing_rate,
    half_saturation,
    palatability,
    (reference_inventory,),
    (palatability,),
    1,
)

end # module
