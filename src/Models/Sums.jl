"""Shared reduction helpers.

This module centralizes small, kernel-friendly reductions that are reused
across model tracer kernels (e.g. NiPiZD and DARWIN).

"""
module Sums

using ...Utils: sum_over, TendencyContext

using ...Library.Mortality: linear_loss, quadratic_loss
using ...Library.Photosynthesis: smith_single_nutrient_growth, geider_two_nutrient_growth
using ...Library.Predation:
    preferential_predation_loss,
    preferential_predation_gain,
    preferential_predation_unassimilated_loss

# -----------------------------------------------------------------------------
# Predation / grazing matrix reductions
# -----------------------------------------------------------------------------

"""Sum of unassimilated grazing loss terms across all consumer-by-prey pairs."""
@inline function grazing_unassimilated_loss_sum(tendency::TendencyContext)
    parameters = tendency.parameters
    tracers = tendency.tracers
    args = tendency.args

    FT = eltype(parameters.maximum_predation_rate)
    z = zero(FT)

    ints = parameters.interactions
    pal = ints.palatability
    assim = ints.assimilation
    consumer_global = ints.consumer_global
    prey_global = ints.prey_global

    sum_over(eachindex(consumer_global), z) do ic
        predator_idx = consumer_global[ic]
        predator = tracers.plankton(args, predator_idx)
        gmax = parameters.maximum_predation_rate[predator_idx]
        K = parameters.holling_half_saturation[predator_idx]

        sum_over(eachindex(prey_global), z) do ip
            prey_idx = prey_global[ip]
            prey = tracers.plankton(args, prey_idx)
            beta = assim[ic, ip]
            phi = pal[ic, ip]
            preferential_predation_unassimilated_loss(prey, predator, beta, gmax, K, phi)
        end
    end
end

"""Sum of grazing loss terms for a single prey (given global index)."""
@inline function grazing_loss_sum(
    tendency::TendencyContext,
    prey,
    prey_idx::Int,
)
    parameters = tendency.parameters
    tracers = tendency.tracers
    args = tendency.args

    FT = eltype(parameters.maximum_predation_rate)
    z = zero(FT)

    ints = parameters.interactions
    ip = @inbounds ints.global_to_prey[prey_idx]
    ip == 0 && return z

    pal = ints.palatability
    consumer_global = ints.consumer_global

    sum_over(eachindex(consumer_global), z) do ic
        predator_idx = consumer_global[ic]
        predator = tracers.plankton(args, predator_idx)
        gmax = parameters.maximum_predation_rate[predator_idx]
        K = parameters.holling_half_saturation[predator_idx]
        phi = pal[ic, ip]
        preferential_predation_loss(prey, predator, gmax, K, phi)
    end
end

"""Sum of grazing gain terms for a single predator (given global index)."""
@inline function grazing_gain_sum(
    tendency::TendencyContext,
    predator,
    predator_idx::Int,
)
    parameters = tendency.parameters
    tracers = tendency.tracers
    args = tendency.args

    FT = eltype(parameters.maximum_predation_rate)
    z = zero(FT)

    ints = parameters.interactions
    ic = @inbounds ints.global_to_consumer[predator_idx]
    ic == 0 && return z

    pal = ints.palatability
    assim = ints.assimilation
    prey_global = ints.prey_global

    gmax = parameters.maximum_predation_rate[predator_idx]
    K = parameters.holling_half_saturation[predator_idx]

    sum_over(eachindex(prey_global), z) do ip
        prey_idx = prey_global[ip]
        prey = tracers.plankton(args, prey_idx)
        beta = assim[ic, ip]
        phi = pal[ic, ip]
        preferential_predation_gain(prey, predator, beta, gmax, K, phi)
    end
end

# -----------------------------------------------------------------------------
# Mortality reductions (plankton-class sums)
# -----------------------------------------------------------------------------

@inline function linear_mortality_sum(n_plankton::Int, vals, linear_mortality)
    sum_over(n_plankton, zero(eltype(linear_mortality))) do i
        P = vals.plankton(i)
        linear_loss(P, linear_mortality[i])
    end
end

@inline function quadratic_mortality_sum(n_plankton::Int, vals, quadratic_mortality)
    sum_over(n_plankton, zero(eltype(quadratic_mortality))) do i
        P = vals.plankton(i)
        quadratic_loss(P, quadratic_mortality[i])
    end
end

@inline function mortality_loss_sum(n_plankton::Int, vals, linear_mortality, quadratic_mortality)
    FT = promote_type(eltype(linear_mortality), eltype(quadratic_mortality))
    z = zero(FT)
    sum_over(n_plankton, z) do i
        P = vals.plankton(i)
        linear_loss(P, linear_mortality[i]) + quadratic_loss(P, quadratic_mortality[i])
    end
end

# -----------------------------------------------------------------------------
# Uptake reductions (growth summed over plankton classes)
# -----------------------------------------------------------------------------

"""Smith-style single-nutrient uptake summed over all plankton classes."""
@inline function smith_uptake_sum(
    n_plankton::Int,
    vals,
    N,
    PAR,
    maximum_growth_rate,
    nutrient_half_saturation,
    alpha,
)
    sum_over(n_plankton, zero(N)) do i
        P = vals.plankton(i)
        smith_single_nutrient_growth(
            N,
            P,
            PAR,
            maximum_growth_rate[i],
            nutrient_half_saturation[i],
            alpha[i],
        )
    end
end

"""Geider-style two-nutrient uptake summed over all plankton classes."""
@inline function geider_two_nutrient_uptake_sum(
    n_plankton::Int,
    vals,
    DIN,
    PO4,
    PAR,
    maximum_growth_rate,
    half_saturation_DIN,
    half_saturation_PO4,
    photosynthetic_slope,
    chlorophyll_to_carbon_ratio,
)
    sum_over(n_plankton, zero(DIN)) do i
        P = vals.plankton(i)
        geider_two_nutrient_growth(
            DIN,
            PO4,
            P,
            PAR,
            maximum_growth_rate[i],
            half_saturation_DIN[i],
            half_saturation_PO4[i],
            photosynthetic_slope[i],
            chlorophyll_to_carbon_ratio[i],
        )
    end
end

end # module Sums
