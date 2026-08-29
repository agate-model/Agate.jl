import ..Parameters: derive_default, _derive_parameter_default

using ..Library.Allometry:
    palatability_matrix_allometric_axes, assimilation_efficiency_matrix_binary_axes

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

"""Derive consumer-by-prey palatability from allometric trait vectors."""
struct PalatabilityAllometric end

"""Derive consumer-by-prey assimilation from binary efficiency traits."""
struct AssimilationBinary end

function _size_class_indices(
    layout::ModelLayout, labels::Tuple, parameter_name::Symbol, axis_name::Symbol
)
    return map(labels) do label
        index = findfirst(==(label), layout.size_classes)
        isnothing(index) && throw(ArgumentError(
            "derived default :$parameter_name $axis_name SizeClass :$label is not a " *
            "realized plankton SizeClass",
        ))
        index
    end
end

@inline function _derive_palatability(layout::ModelLayout, params, consumers, prey)
    T = layout.scalar_type
    return palatability_matrix_allometric_axes(
        T,
        layout.size_class_diameters;
        optimum_predator_prey_ratio=_require_scalar_vector(
            T, params.optimum_predator_prey_ratio, :optimum_predator_prey_ratio
        ),
        specificity=_require_scalar_vector(T, params.specificity, :specificity),
        protection=_require_scalar_vector(T, params.protection, :protection),
        consumer_indices=consumers,
        prey_indices=prey,
    )
end

@inline function _derive_assimilation(layout::ModelLayout, params, consumers, prey)
    T = layout.scalar_type
    return assimilation_efficiency_matrix_binary_axes(
        T;
        assimilation_efficiency=_require_scalar_vector(
            T, params.assimilation_efficiency, :assimilation_efficiency
        ),
        consumer_indices=consumers,
        prey_indices=prey,
    )
end

@inline function derive_default(
    ::PalatabilityAllometric, ::Any, layout::ModelLayout, params::NamedTuple
)
    return _derive_palatability(
        layout, params, layout.consumer_indices, layout.prey_indices
    )
end

@inline function derive_default(
    ::AssimilationBinary, ::Any, layout::ModelLayout, params::NamedTuple
)
    return _derive_assimilation(
        layout, params, layout.consumer_indices, layout.prey_indices
    )
end

@inline function _derive_parameter_default(
    ::PalatabilityAllometric,
    ::Any,
    layout::ModelLayout,
    parameter,
    params::NamedTuple,
)
    consumer_labels, prey_labels = parameter.storage_labels
    return _derive_palatability(
        layout,
        params,
        _size_class_indices(layout, consumer_labels, parameter.name, :consumer),
        _size_class_indices(layout, prey_labels, parameter.name, :prey),
    )
end

@inline function _derive_parameter_default(
    ::AssimilationBinary,
    ::Any,
    layout::ModelLayout,
    parameter,
    params::NamedTuple,
)
    consumer_labels, resource_labels = parameter.storage_labels
    return _derive_assimilation(
        layout,
        params,
        _size_class_indices(layout, consumer_labels, parameter.name, :consumer),
        eachindex(resource_labels),
    )
end
