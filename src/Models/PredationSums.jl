"""Shared predation reduction helpers.

These functions implement the consumer-by-prey rectangular interaction matrix
pattern using the canonical `parameters.interactions` container.

They are shared by multiple model tracer kernels (e.g. NiPiZD and DARWIN) to
avoid duplicated indexing logic.

GPU notes
---------
All functions operate on positional tracer arguments via a `Tracers` accessor.
No runtime `Symbol` lookup is performed in kernel-callable code.
"""
module PredationSums

using ...Utils: sum_over, TendencyContext
using ...Library.Predation:
    preferential_predation_loss,
    preferential_predation_gain,
    preferential_predation_unassimilated_loss

"""Sum of assimilation loss terms across all consumer-by-prey pairs.

The accumulator is always the model float type `FT`, which is inferred from
the interaction parameters. This keeps kernel call sites simple:

    loss = _grazing_assimilation_loss_sum(tendency)

No generic reduction seed is required because Agate always runs with a single
adapted floating-point type.
"""
@inline function _grazing_assimilation_loss_sum(tendency::TendencyContext)
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
@inline function _grazing_loss_sum(
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
@inline function _grazing_gain_sum(
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

end # module PredationSums
