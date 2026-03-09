"""Abstract supertype for strategies that derive interaction matrices.

Derivers run during factory construction on the host and must return concrete
matrices that can be adapted to the target architecture with `Adapt`.
"""
abstract type AbstractMatrixDeriver end

"""Describe how to build a single derived matrix.

Fields
------
- `deriver`: strategy object that computes the matrix.
- `deps`: parameter keys that trigger recomputation when explicitly overridden.
"""
struct MatrixDefinition{D<:AbstractMatrixDeriver}
    deriver::D
    deps::Tuple{Vararg{Symbol}}
end

function MatrixDefinition(
    deriver::D; deps=derivation_deps(deriver)
) where {D<:AbstractMatrixDeriver}
    deps isa Tuple{Vararg{Symbol}} ||
        throw(ArgumentError("deps must be a tuple of Symbols (got $(typeof(deps)))."))
    return MatrixDefinition{D}(deriver, deps)
end

"""Return the dependency keys for a derivation strategy.

Factories use these keys to decide whether a derived matrix should be
recomputed when inputs are explicitly overridden.
"""
derivation_deps(::AbstractMatrixDeriver) = ()

"""Compute a derived matrix for `deriver`.

Concrete derivers must implement
`derive_matrix(::Deriver, factory, context, params)`.
"""
function derive_matrix(
    ::AbstractMatrixDeriver, ::AbstractBGCFactory, ::CommunityContext, ::NamedTuple
)
    throw(ArgumentError("derive_matrix is not implemented for this deriver"))
end

"""Return derivation specs for matrices exposed by `factory`.

Some models expose low-level interaction traits (for example,
predator:prey size-ratio preferences) and derive higher-level interaction
matrices such as palatability or assimilation from those traits.

Derived matrices are computed during construction on the host when they are
missing. Existing values are recomputed only when both of the following are
true:

- the derived matrix itself was not explicitly overridden, and
- at least one declared dependency was explicitly overridden.

This lets advanced users override traits while keeping user-facing interaction
overrides data-only, as explicit rectangular matrices.

Advanced usage
--------------

Factories declare derivable matrices by returning a `NamedTuple` of
`MatrixDefinition` objects from `matrix_definitions(factory)`. Each entry maps a
parameter key to a derivation strategy and its dependency keys.

To change the derivation algorithm for a matrix, swap the deriver in
`matrix_definitions`. For example:

```julia
function matrix_definitions(::MyFactory)
    return (;
        assimilation_matrix = MatrixDefinition(AssimilationNewFunction()),
    )
end
```

`resolve_derived_matrices` enforces that every derived matrix is a concrete
rectangular matrix with canonical axes and `eltype == FT`, where `FT` is the
grid floating-point type. No implicit casting is performed.
"""
matrix_definitions(::AbstractBGCFactory) = (;)

"""Validate the shape and element type of a derived matrix result."""
function _validate_derived_matrix_result(
    factory::AbstractBGCFactory, context::CommunityContext, key::Symbol, value
)
    spec = parameter_spec(factory, key)
    spec !== nothing || throw(
        ArgumentError(
            "derived matrix '$key' is missing a ParameterSpec in parameter_directory(::$(typeof(factory))).",
        ),
    )

    spec.shape === :matrix || throw(
        ArgumentError(
            "derived matrix '$key' must target a matrix ParameterSpec; got shape $(spec.shape).",
        ),
    )

    value isa AbstractMatrix || throw(
        ArgumentError("derived matrix '$key' must be a matrix; got $(typeof(value)).")
    )

    FT = context.FT
    eltype(value) === FT || throw(
        ArgumentError(
            "derived matrix '$key' must have eltype $(FT); got eltype $(eltype(value)) (type: $(typeof(value))). No implicit casting is performed.",
        ),
    )

    if spec.axes === nothing
        n = context.n_total
        size(value) == (n, n) || throw(
            ArgumentError(
                "derived matrix '$key' must be a $(n)x$(n) matrix; got size $(size(value)).",
            ),
        )
        return nothing
    end

    row_axis, col_axis = spec.axes
    nr = length(axis_indices(context, row_axis))
    nc = length(axis_indices(context, col_axis))

    size(value) == (nr, nc) || throw(
        ArgumentError(
            "derived matrix '$key' expected rectangular matrix of size ($(nr), $(nc)) for axes $(spec.axes); got size $(size(value)).",
        ),
    )

    return nothing
end

"""Resolve derived matrices into `params`.

`explicit_override_keys` must list keys explicitly provided by the user, either
through `parameters` or `interaction_overrides`.

For each derived matrix key `K`, resolution follows this order:

1. If `K` was explicitly overridden, keep the provided value.
2. Otherwise, if `K` is missing, compute it.
3. Otherwise, if any declared dependency was explicitly overridden,
   recompute it.
4. Otherwise, keep the existing value.
"""
function resolve_derived_matrices(
    factory::AbstractBGCFactory,
    context::CommunityContext,
    params::NamedTuple,
    explicit_override_keys::Tuple{Vararg{Symbol}},
)
    derived_specs = matrix_definitions(factory)
    isempty(derived_specs) && return params

    override_set = Set(explicit_override_keys)
    resolved = params

    for key in propertynames(derived_specs)
        spec = getproperty(derived_specs, key)
        spec isa MatrixDefinition || throw(
            ArgumentError(
                "matrix_definitions(::$(typeof(factory))) must return MatrixDefinition entries; got $(typeof(spec)) for key :$key.",
            ),
        )

        # Rule 1: explicit override wins.
        if (key in override_set) && hasproperty(resolved, key)
            continue
        end

        needs_compute = !hasproperty(resolved, key)
        needs_recompute = !isempty(spec.deps) && any(d -> d in override_set, spec.deps)

        if needs_compute || needs_recompute
            value = derive_matrix(spec.deriver, factory, context, resolved)
            _validate_derived_matrix_result(factory, context, key, value)
            resolved = merge(resolved, NamedTuple{(key,)}((value,)))
        end
    end

    return resolved
end

# -----------------------------------------------------------------------------
# Built-in derivers
# -----------------------------------------------------------------------------

using ..Library.Allometry:
    palatability_matrix_allometric_axes, assimilation_efficiency_matrix_binary_axes

"""Return `v` if it has element type `FT`, otherwise throw an `ArgumentError`."""
@inline function _require_FT_vector(
    ::Type{FT}, v::AbstractVector, name::Symbol
) where {FT<:AbstractFloat}
    eltype(v) === FT && return v
    throw(
        ArgumentError(
            "expected trait vector :$name to have eltype $(FT); got eltype $(eltype(v)) (type: $(typeof(v))). " *
            "Traits must be provided as Vector{$(FT)} (the grid float type). " *
            "No implicit casting is performed; define traits in a Variant/Factory default.",
        ),
    )
end

"""Derive `palatability_matrix[consumer, prey]` from allometric trait vectors."""
struct PalatabilityAllometric <: AbstractMatrixDeriver end

"""Derive `assimilation_matrix[consumer, prey]` from binary efficiency traits."""
struct AssimilationBinary <: AbstractMatrixDeriver end

# Dependency keys are the trait vectors used by each derivation.
function derivation_deps(::PalatabilityAllometric)
    return (:optimum_predator_prey_ratio, :specificity, :protection)
end
derivation_deps(::AssimilationBinary) = (:assimilation_efficiency,)

@inline function derive_matrix(
    ::PalatabilityAllometric,
    ::AbstractBGCFactory,
    context::CommunityContext,
    params::NamedTuple,
)
    FT = context.FT
    return palatability_matrix_allometric_axes(
        FT,
        context.diameters;
        optimum_predator_prey_ratio=_require_FT_vector(
            FT, params.optimum_predator_prey_ratio, :optimum_predator_prey_ratio
        ),
        specificity=_require_FT_vector(FT, params.specificity, :specificity),
        protection=_require_FT_vector(FT, params.protection, :protection),
        consumer_indices=context.consumer_indices,
        prey_indices=context.prey_indices,
    )
end

@inline function derive_matrix(
    ::AssimilationBinary,
    ::AbstractBGCFactory,
    context::CommunityContext,
    params::NamedTuple,
)
    FT = context.FT
    return assimilation_efficiency_matrix_binary_axes(
        FT;
        assimilation_efficiency=_require_FT_vector(
            FT, params.assimilation_efficiency, :assimilation_efficiency
        ),
        consumer_indices=context.consumer_indices,
        prey_indices=context.prey_indices,
    )
end
