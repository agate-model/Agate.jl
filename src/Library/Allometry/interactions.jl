# -----------------------------------------------------------------------------
# Palatability + assimilation matrix utilities
# -----------------------------------------------------------------------------

"""
    PalatabilityPreyParameters(diameter, protection)

Prey traits used by allometric palatability kernels.

# Fields
- `diameter`: prey equivalent spherical diameter.
- `protection`: dimensionless prey protection factor, where `0` means no
  protection and `1` means complete protection in the protected palatability
  kernel.
"""
struct PalatabilityPreyParameters{T<:Real}
    diameter::T
    protection::T
end

"""
    PalatabilityPredatorParameters(diameter, optimum_predator_prey_ratio, specificity)

Predator traits used by allometric palatability kernels.

# Fields
- `diameter`: predator equivalent spherical diameter.
- `optimum_predator_prey_ratio`: preferred predator:prey diameter ratio.
- `specificity`: sharpness of the unimodal prey-size preference.
"""
struct PalatabilityPredatorParameters{T<:Real}
    diameter::T
    optimum_predator_prey_ratio::T
    specificity::T
end

"""
    allometric_scaling_power(a, b, diameter)

Evaluate a power-law allometric scaling against spherical cell volume.

!!! formulation
    ```math
    f(d) = a V(d)^b,
    \\qquad
    V(d) = \\frac{4}{3}\\pi\\left(\\frac{d}{2}\\right)^3
    ```

    where ``d`` is equivalent spherical diameter, ``V`` is spherical cell volume,
    ``a`` is the prefactor, and ``b`` is the exponent.

# Arguments
- `a`: scale/prefactor parameter.
- `b`: exponent parameter.
- `diameter`: cell equivalent spherical diameter.
"""
@inline function allometric_scaling_power(a::T, b::T, diameter::T) where {T<:Real}
    r = diameter / T(2)
    volume = (T(4) / T(3)) * T(π) * r^T(3)
    return a * volume^b
end

"""
    allometric_palatability_unimodal(prey, predator)

Compute unimodal allometric palatability from predator and prey diameters.

!!! formulation
    ```math
    \\eta = \\left[1 + \\left(\\frac{d_{pred}}{d_{prey}} - \\rho^*\\right)^2\\right]^{-\\sigma}
    ```

    where ``d_{pred}`` and ``d_{prey}`` are predator and prey diameters,
    ``\\rho^*`` is the optimum predator:prey diameter ratio, and ``\\sigma`` is
    the specificity parameter.

# Arguments
- `prey`: `PalatabilityPreyParameters(diameter, protection)`.
- `predator`: `PalatabilityPredatorParameters(diameter, optimum_predator_prey_ratio, specificity)`.
"""
@inline function allometric_palatability_unimodal(
    prey::PalatabilityPreyParameters{T}, predator::PalatabilityPredatorParameters{T}
) where {T<:Real}
    ratio = predator.diameter / prey.diameter
    width = one(T) + (ratio - predator.optimum_predator_prey_ratio)^T(2)
    return one(T) / width^predator.specificity
end

"""
    allometric_palatability_unimodal_protection(prey, predator)

Compute unimodal allometric palatability with multiplicative prey protection.

!!! formulation
    ```math
    \\eta = (1 - p)\\left[1 + \\left(\\frac{d_{pred}}{d_{prey}} - \\rho^*\\right)^2\\right]^{-\\sigma}
    ```

    where ``p`` is `prey.protection`. `p = 0` leaves the allometric
    palatability unchanged, while larger values reduce palatability.
"""
function allometric_palatability_unimodal_protection(
    prey::PalatabilityPreyParameters{T}, predator::PalatabilityPredatorParameters{T}
) where {T<:Real}
    base = allometric_palatability_unimodal(prey, predator)
    return base * (one(T) - prey.protection)
end

"""
    palatability_matrix_allometric_axes(T, diameters; optimum_predator_prey_ratio,
                                        specificity, protection,
                                        consumer_indices, prey_indices,
                                        palatability_fn=allometric_palatability_unimodal_protection)

Build a consumer-by-prey palatability matrix from allometric traits.

!!! formulation
    ```math
    M_{ij} = \\eta(prey_j, consumer_i)
    ```

    where ``\\eta`` is `palatability_fn`. Only rows from `consumer_indices` and
    columns from `prey_indices` are materialized, so the returned matrix has size
    `length(consumer_indices) × length(prey_indices)`.

# Arguments
- `T`: scalar type used for the output matrix and trait values.
- `diameters`: full diameter vector indexed by model class.
- `optimum_predator_prey_ratio`: full vector of predator optimum ratios.
- `specificity`: full vector of predator specificity parameters.
- `protection`: full vector of prey protection parameters.
- `consumer_indices`: source indices for matrix rows.
- `prey_indices`: source indices for matrix columns.
- `palatability_fn`: callable mapping prey and predator trait structs to a scalar.
"""
function palatability_matrix_allometric_axes(
    ::Type{T},
    diameters::AbstractVector{T};
    optimum_predator_prey_ratio::AbstractVector{T},
    specificity::AbstractVector{T},
    protection::AbstractVector{T},
    consumer_indices::AbstractVector{<:Integer},
    prey_indices::AbstractVector{<:Integer},
    palatability_fn=allometric_palatability_unimodal_protection,
) where {T<:Real}
    nr = length(consumer_indices)
    nc = length(prey_indices)
    M = zeros(T, nr, nc)

    @inbounds for (ii, pred) in pairs(consumer_indices)
        predator = PalatabilityPredatorParameters{T}(
            diameters[pred], optimum_predator_prey_ratio[pred], specificity[pred]
        )
        for (jj, prey) in pairs(prey_indices)
            prey_params = PalatabilityPreyParameters{T}(diameters[prey], protection[prey])
            M[ii, jj] = palatability_fn(prey_params, predator)
        end
    end

    return M
end

"""
    assimilation_efficiency_matrix_binary_axes(T; assimilation_efficiency,
                                               consumer_indices, prey_indices)

Build a consumer-by-prey assimilation-efficiency matrix.

!!! formulation
    ```math
    B_{ij} = \\beta_i
    ```

    where ``\\beta_i`` is the assimilation efficiency of consumer `i`. Only rows
    from `consumer_indices` and columns from `prey_indices` are materialized, so
    the returned matrix has size `length(consumer_indices) × length(prey_indices)`.

# Arguments
- `T`: scalar type used for the output matrix.
- `assimilation_efficiency`: full vector of consumer assimilation efficiencies.
- `consumer_indices`: source indices for matrix rows.
- `prey_indices`: source indices for matrix columns.
"""
function assimilation_efficiency_matrix_binary_axes(
    ::Type{T};
    assimilation_efficiency::AbstractVector{T},
    consumer_indices::AbstractVector{<:Integer},
    prey_indices::AbstractVector{<:Integer},
) where {T<:Real}
    nr = length(consumer_indices)
    nc = length(prey_indices)
    M = zeros(T, nr, nc)

    @inbounds for (ii, pred) in pairs(consumer_indices)
        β = assimilation_efficiency[pred]
        for jj in 1:nc
            M[ii, jj] = β
        end
    end

    return M
end
