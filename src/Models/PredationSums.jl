"""Shared predation reduction helpers.

These functions implement the consumer-by-prey rectangular interaction matrix
pattern using the canonical `parameters.interactions` container.

They are shared by multiple model tracer kernels (e.g. NiPiZD and DARWIN) to
avoid duplicated indexing logic.
"""
module PredationSums

using ...Utils: sum_over
using ...Library.Predation:
    PreferentialPredationLoss,
    PreferentialPredationGain,
    PreferentialPredationAssimilationLoss

"""Sum of assimilation loss terms across all consumer-by-prey pairs.

Returns a value with the same type as `init`.
"""
@inline function _grazing_assimilation_loss_sum(p, state, plankton_syms, init)
    ints = p.interactions
    pal = ints.palatability
    assim = ints.assimilation
    consumer_global = ints.consumer_global
    prey_global = ints.prey_global

    sum_over(eachindex(consumer_global), init) do ic
        predator_idx = consumer_global[ic]
        predator = getproperty(state, plankton_syms[predator_idx])
        gmax = p.maximum_predation_rate[predator_idx]
        K = p.holling_half_saturation[predator_idx]

        sum_over(eachindex(prey_global), zero(init)) do ip
            prey_idx = prey_global[ip]
            prey = getproperty(state, plankton_syms[prey_idx])
            beta = assim[ic, ip]
            phi = pal[ic, ip]
            PreferentialPredationAssimilationLoss(beta, gmax, K, phi)(prey, predator)
        end
    end
end

"""Sum of grazing loss terms for a single prey (given global index)."""
@inline function _grazing_loss_sum(p, state, plankton_syms, prey, prey_idx::Int, init)
    ints = p.interactions
    ip = @inbounds ints.global_to_prey[prey_idx]
    ip == 0 && return init

    pal = ints.palatability
    consumer_global = ints.consumer_global

    sum_over(eachindex(consumer_global), init) do ic
        predator_idx = consumer_global[ic]
        predator = getproperty(state, plankton_syms[predator_idx])
        gmax = p.maximum_predation_rate[predator_idx]
        K = p.holling_half_saturation[predator_idx]
        phi = pal[ic, ip]
        PreferentialPredationLoss(gmax, K, phi)(prey, predator)
    end
end

"""Sum of grazing gain terms for a single predator (given global index)."""
@inline function _grazing_gain_sum(p, state, plankton_syms, predator, predator_idx::Int, init)
    ints = p.interactions
    ic = @inbounds ints.global_to_consumer[predator_idx]
    ic == 0 && return init

    pal = ints.palatability
    assim = ints.assimilation
    prey_global = ints.prey_global

    gmax = p.maximum_predation_rate[predator_idx]
    K = p.holling_half_saturation[predator_idx]

    sum_over(eachindex(prey_global), init) do ip
        prey_idx = prey_global[ip]
        prey = getproperty(state, plankton_syms[prey_idx])
        beta = assim[ic, ip]
        phi = pal[ic, ip]
        PreferentialPredationGain(beta, gmax, K, phi)(prey, predator)
    end
end

end # module PredationSums
