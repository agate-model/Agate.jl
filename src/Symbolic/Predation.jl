module Predation

using ...Equations: sum_over

export grazing_loss, grazing_gain, grazing_assimilation_loss

@inline _vec(PV, key::Symbol, idx::Int) = getproperty(PV, key)[idx]
@inline _mat(PV, key::Symbol, j::Int, i::Int) = getproperty(PV, key)[j, i]

"""
    grazing_loss(PV, prey_sym, prey_idx, plankton_syms; ...)

Preferential grazing loss of a prey size-class to all potential predators.
"""
function grazing_loss(
    PV,
    prey_sym::Symbol,
    prey_idx::Int,
    plankton_syms::AbstractVector{Symbol};
    rate::Symbol=:maximum_predation_rate,
    half::Symbol=:holling_half_saturation,
    palat::Symbol=:palatability_matrix,
)
    return sum_over(plankton_syms) do predator_sym, predator_idx
        gmax = _vec(PV, rate, predator_idx)
        K = _vec(PV, half, predator_idx)
        ϕ = _mat(PV, palat, predator_idx, prey_idx)
        return gmax * ϕ * (prey_sym / (K + prey_sym)) * predator_sym
    end
end

"""
    grazing_gain(PV, predator_sym, predator_idx, plankton_syms; ...)

Preferential grazing gain of a predator size-class from all potential prey.
"""
function grazing_gain(
    PV,
    predator_sym::Symbol,
    predator_idx::Int,
    plankton_syms::AbstractVector{Symbol};
    assim::Symbol=:assimilation_matrix,
    rate::Symbol=:maximum_predation_rate,
    half::Symbol=:holling_half_saturation,
    palat::Symbol=:palatability_matrix,
)
    return sum_over(plankton_syms) do prey_sym, prey_idx
        β = _mat(PV, assim, predator_idx, prey_idx)
        gmax = _vec(PV, rate, predator_idx)
        K = _vec(PV, half, predator_idx)
        ϕ = _mat(PV, palat, predator_idx, prey_idx)
        return β * gmax * ϕ * (prey_sym / (K + prey_sym)) * predator_sym
    end
end

"""
    grazing_assimilation_loss(PV, plankton_syms; ...)

Sum unassimilated grazing losses across all predator–prey pairs.
"""
function grazing_assimilation_loss(
    PV,
    plankton_syms::AbstractVector{Symbol};
    assim::Symbol=:assimilation_matrix,
    rate::Symbol=:maximum_predation_rate,
    half::Symbol=:holling_half_saturation,
    palat::Symbol=:palatability_matrix,
)
    return sum_over(plankton_syms) do predator_sym, predator_idx
        sum_over(plankton_syms) do prey_sym, prey_idx
            β = _mat(PV, assim, predator_idx, prey_idx)
            gmax = _vec(PV, rate, predator_idx)
            K = _vec(PV, half, predator_idx)
            ϕ = _mat(PV, palat, predator_idx, prey_idx)
            return (1 - β) * gmax * ϕ * (prey_sym / (K + prey_sym)) * predator_sym
        end
    end
end

end # module
