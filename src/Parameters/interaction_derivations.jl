using ..Components: ModelLayout, diameter_metadata

using ..Library.Allometry: palatability_matrix_allometric_axes

"""Return `v` when it uses the construction scalar type, otherwise throw an `ArgumentError`."""
@inline function _require_scalar_vector(
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

"""Derive consumer-by-prey palatability from size traits, with optional prey protection."""
struct AllometricPalatability end

"""Broadcast one consumer-indexed value across every resource for that consumer.

The single declared dependency must be a vector over realized plankton SizeClasses. The
result is a consumer-by-resource matrix whose rows repeat the corresponding consumer value.
"""
struct BroadcastConsumerAcrossResources end

function _plankton_entity_indices(
    layout::ModelLayout, labels::Tuple, parameter_name::Symbol, axis_name::Symbol
)
    return map(labels) do label
        index = findfirst(==(label), layout.size_classes)
        isnothing(index) && throw(
            ArgumentError(
                "derived default :$parameter_name $axis_name entity :$label is not a " *
                "realized plankton entity",
            ),
        )
        return index
    end
end

function _require_palatability_diameters(layout::ModelLayout, consumers, prey)
    seen = Set{Int}()
    for index in (consumers..., prey...)
        index in seen && continue
        push!(seen, index)
        diameter_metadata(layout.size_class_diameters[index]) === nothing || continue
        entity = layout.size_classes[index]
        throw(
            ArgumentError(
                "AllometricPalatability requires diameter metadata for SizeClass :$entity; provide an explicit size structure or an explicit palatability matrix",
            ),
        )
    end
    return nothing
end

# Prepare the fixed numeric context once; both construction and trait evaluation
# invoke the same allometric matrix kernel.
function _prepare_palatability(::Type{T}, diameters, consumers, prey) where {T<:Real}
    return function (params)
        protection = if hasproperty(params, :protection)
            _require_scalar_vector(T, params.protection, :protection)
        else
            zeros(T, length(diameters))
        end
        return palatability_matrix_allometric_axes(
            T,
            diameters;
            optimum_predator_prey_ratio=_require_scalar_vector(
                T, params.optimum_predator_prey_ratio, :optimum_predator_prey_ratio
            ),
            specificity=_require_scalar_vector(T, params.specificity, :specificity),
            protection,
            consumer_indices=consumers,
            prey_indices=prey,
        )
    end
end

function _prepare_parameter_default(
    ::AllometricPalatability, ::Any, layout::ModelLayout, parameter
)
    consumer_labels, prey_labels = parameter.storage_labels
    consumers = _plankton_entity_indices(layout, consumer_labels, parameter.name, :consumer)
    prey = _plankton_entity_indices(layout, prey_labels, parameter.name, :prey)
    _require_palatability_diameters(layout, consumers, prey)
    return _prepare_palatability(
        layout.scalar_type, layout.size_class_diameters, consumers, prey
    )
end

function _prepare_consumer_broadcast(::Type{T}, consumers, nresources) where {T<:Real}
    return function (params)
        length(params) == 1 || throw(
            ArgumentError(
                "BroadcastConsumerAcrossResources requires exactly one consumer-indexed dependency",
            ),
        )
        trait = _require_scalar_vector(T, first(values(params)), first(keys(params)))
        return [trait[i] for i in consumers, _ in 1:nresources]
    end
end

function _prepare_parameter_default(
    ::BroadcastConsumerAcrossResources, ::Any, layout::ModelLayout, parameter
)
    consumer_labels, resource_labels = parameter.storage_labels
    consumers = _plankton_entity_indices(layout, consumer_labels, parameter.name, :consumer)
    return _prepare_consumer_broadcast(
        layout.scalar_type, consumers, length(resource_labels)
    )
end
