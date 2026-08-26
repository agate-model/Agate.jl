import ..Parameters: derive_default

using ..Library.Allometry:
    palatability_matrix_allometric_axes, assimilation_efficiency_matrix_binary_axes

"""Return `v` when it uses the construction scalar type, otherwise throw an `ArgumentError`."""
@inline function require_scalar_vector(
    ::Type{T}, v::AbstractVector, name::Symbol
) where {T<:Real}
    eltype(v) === T && return v
    throw(
        ArgumentError(
            "expected trait vector :$name to have eltype $(T); got eltype $(eltype(v)) (type: $(typeof(v))). " *
            "Traits must use the construction scalar type. No implicit casting is performed.",
        ),
    )
end

"""Derive consumer-by-prey palatability from allometric trait vectors."""
struct PalatabilityAllometric end

"""Derive consumer-by-prey assimilation from binary efficiency traits."""
struct AssimilationBinary end

@inline function derive_default(
    ::PalatabilityAllometric,
    ::Any,
    layout::ModelLayout,
    params::NamedTuple,
)
    T = layout.scalar_type
    return palatability_matrix_allometric_axes(
        T,
        layout.diameters;
        optimum_predator_prey_ratio=require_scalar_vector(
            T, params.optimum_predator_prey_ratio, :optimum_predator_prey_ratio
        ),
        specificity=require_scalar_vector(T, params.specificity, :specificity),
        protection=require_scalar_vector(T, params.protection, :protection),
        consumer_indices=layout.consumer_indices,
        prey_indices=layout.prey_indices,
    )
end

@inline function derive_default(
    ::AssimilationBinary,
    ::Any,
    layout::ModelLayout,
    params::NamedTuple,
)
    T = layout.scalar_type
    return assimilation_efficiency_matrix_binary_axes(
        T;
        assimilation_efficiency=require_scalar_vector(
            T, params.assimilation_efficiency, :assimilation_efficiency
        ),
        consumer_indices=layout.consumer_indices,
        prey_indices=layout.prey_indices,
    )
end
