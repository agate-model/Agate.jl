"""Shared reduction helpers.

This module centralizes small, kernel-friendly reductions that are reused
across model tracer kernels (e.g. NiPiZD and DARWIN).

The helpers here are **domain-specific** (Agate model internals) and are
intended to support a hot-path tracer DSL with predictable types.
"""
module Sums

using ...Library.Mortality: linear_loss, quadratic_loss
using ...Library.Photosynthesis: smith_single_nutrient_growth, geider_two_nutrient_growth
using ...Library.Predation:
    preferential_predation_loss,
    preferential_predation_gain,
    preferential_predation_unassimilated_loss

# -----------------------------------------------------------------------------
# Local reduction helper
# -----------------------------------------------------------------------------

"""Sum `f(i)` for `i` in `itr` with accumulator type `FT`.

Designed for Julia's `do`-block syntax, e.g.

```julia
sum_over(n, FT) do i
    ...
end
```

`itr` may be an integer `n` (interpreted as `1:n`), a range, or `eachindex(x)`.
"""
@inline function sum_over(f, itr, ::Type{FT}) where {FT}
    acc = zero(FT)
    @inbounds for i in itr
        acc += f(i)
    end
    return acc
end

@inline sum_over(f, n::Int, ::Type{FT}) where {FT} = sum_over(f, Base.OneTo(n), FT)

# -----------------------------------------------------------------------------
# Predation / grazing matrix reductions
# -----------------------------------------------------------------------------

"""Sum of unassimilated grazing loss terms across all consumer-by-prey pairs."""
function grazing_unassimilated_loss_sum(parameters, plankton)
    ints = parameters.interactions
    consumer_global = ints.consumer_global
    prey_global = ints.prey_global

    FT = eltype(parameters.maximum_predation_rate)

    n_cons = length(consumer_global)
    n_prey = length(prey_global)

    if n_cons == 0 || n_prey == 0
        return zero(FT)
    end

    pal = ints.palatability
    assim = ints.assimilation

    @inbounds sum_over(eachindex(consumer_global), FT) do ic
        predator_idx = consumer_global[ic]
        predator = plankton(predator_idx)
        gmax = parameters.maximum_predation_rate[predator_idx]
        K = parameters.holling_half_saturation[predator_idx]

        @inbounds sum_over(eachindex(prey_global), FT) do ip
            prey_idx = prey_global[ip]
            prey = plankton(prey_idx)
            beta = assim[ic, ip]
            phi = pal[ic, ip]
            preferential_predation_unassimilated_loss(prey, predator, beta, gmax, K, phi)
        end
    end
end

"""Sum of grazing loss terms for a single prey (given global plankton index)."""
function grazing_loss_sum(
    parameters,
    plankton,
    prey,
    prey_idx::Int,
)
    ints = parameters.interactions
    ip = ints.global_to_prey[prey_idx]

    FT = typeof(prey)
    ip == 0 && return zero(FT)

    pal = ints.palatability
    consumer_global = ints.consumer_global

    @inbounds sum_over(eachindex(consumer_global), FT) do ic
        predator_idx = consumer_global[ic]
        predator = plankton(predator_idx)
        gmax = parameters.maximum_predation_rate[predator_idx]
        K = parameters.holling_half_saturation[predator_idx]
        phi = pal[ic, ip]
        preferential_predation_loss(prey, predator, gmax, K, phi)
    end
end

"""Sum of grazing gain terms for a single predator (given global plankton index)."""
function grazing_gain_sum(
    parameters,
    plankton,
    predator,
    predator_idx::Int,
)
    ints = parameters.interactions
    ic = ints.global_to_consumer[predator_idx]

    FT = typeof(predator)
    ic == 0 && return zero(FT)

    pal = ints.palatability
    assim = ints.assimilation
    prey_global = ints.prey_global

    gmax = parameters.maximum_predation_rate[predator_idx]
    K = parameters.holling_half_saturation[predator_idx]

    @inbounds sum_over(eachindex(prey_global), FT) do ip
        prey_idx = prey_global[ip]
        prey = plankton(prey_idx)
        beta = assim[ic, ip]
        phi = pal[ic, ip]
        preferential_predation_gain(prey, predator, beta, gmax, K, phi)
    end
end

# -----------------------------------------------------------------------------
# Mortality reductions (plankton-class sums)
# -----------------------------------------------------------------------------

function linear_mortality_sum(plankton, linear_mortality)
    n_plankton = length(linear_mortality)
    FT = eltype(linear_mortality)
    @inbounds sum_over(n_plankton, FT) do i
        P = plankton(i)
        linear_loss(P, linear_mortality[i])
    end
end

function quadratic_mortality_sum(plankton, quadratic_mortality)
    n_plankton = length(quadratic_mortality)
    FT = eltype(quadratic_mortality)
    @inbounds sum_over(n_plankton, FT) do i
        P = plankton(i)
        quadratic_loss(P, quadratic_mortality[i])
    end
end

function mortality_loss_sum(plankton, linear_mortality, quadratic_mortality)
    n_plankton = length(linear_mortality)
    FT = eltype(linear_mortality)
    @inbounds sum_over(n_plankton, FT) do i
        P = plankton(i)
        linear_loss(P, linear_mortality[i]) + quadratic_loss(P, quadratic_mortality[i])
    end
end

# -----------------------------------------------------------------------------
# Uptake reductions (growth summed over plankton classes)
# -----------------------------------------------------------------------------

"""Smith-style single-nutrient uptake summed over all plankton classes."""
function smith_uptake_sum(
    plankton,
    N,
    PAR,
    maximum_growth_rate,
    nutrient_half_saturation,
    alpha,
)
    n_plankton = length(maximum_growth_rate)
    FT = eltype(maximum_growth_rate)

    @inbounds sum_over(n_plankton, FT) do i
        P = plankton(i)
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
function geider_two_nutrient_uptake_sum(
    plankton,
    DIN,
    PO4,
    PAR,
    maximum_growth_rate,
    half_saturation_DIN,
    half_saturation_PO4,
    photosynthetic_slope,
    chlorophyll_to_carbon_ratio,
)
    n_plankton = length(maximum_growth_rate)
    FT = eltype(maximum_growth_rate)

    @inbounds sum_over(n_plankton, FT) do i
        P = plankton(i)
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
