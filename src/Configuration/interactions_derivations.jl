"""Utilities for *derived* interaction matrices.

Some models expose low-level interaction traits (for example predator:prey size
ratio preferences) and derive interaction matrices (palatability, assimilation)
from those traits.

Derived matrices are computed during construction on the host when missing. They
are recomputed when:

- a derived matrix is *not* explicitly overridden, and
- at least one of its declared dependencies *is* explicitly overridden.

This allows advanced users to tweak traits while keeping user-facing interaction
overrides data-only (explicit rectangular matrices).

Advanced usage
--------------

Factories may declare derivable matrices via `matrix_definitions(factory)`,
which returns a `NamedTuple` mapping parameter keys to `MatrixDefinition`
instances. Each spec selects a derivation strategy object (`AbstractMatrixDeriver`)
plus a list of dependency keys.

To change the derivation algorithm for a matrix, swap the deriver in
`matrix_definitions`. For example:

```julia
function matrix_definitions(::MyFactory)
    return (;
        assimilation_matrix = MatrixDefinition(AssimilationNewFunction()),
    )
end
```

`resolve_derived_matrices` enforces that derived matrices are concrete
rectangular matrices with canonical axes and have `eltype == FT` (the grid float
type). No implicit casting is performed.
"""

"""Abstract supertype for derived matrix strategies.

Derivers run during construction on the host and must return concrete matrices
that can be moved to the target architecture via `Adapt`.
"""
abstract type AbstractMatrixDeriver end

"""Derivation metadata for one derived matrix.

- `deriver` selects the algorithm.
- `deps` lists parameter keys that trigger recomputation when explicitly overridden.
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

"""Return dependency keys for a derivation strategy.

Factories use dependencies to decide when a derived matrix should be recomputed
if its inputs are explicitly overridden.
"""
derivation_deps(::AbstractMatrixDeriver) = ()

"""Compute a derived matrix.

Derivers must implement `derive_matrix(::Deriver, factory, context, params)`.
"""
function derive_matrix(
    ::AbstractMatrixDeriver, ::AbstractBGCFactory, ::CommunityContext, ::NamedTuple
)
    throw(ArgumentError("derive_matrix is not implemented for this deriver"))
end

"""Return a `NamedTuple` mapping matrix keys to `MatrixDefinition` entries.

Factories override this to declare which matrices can be derived from other
parameters.
"""
matrix_definitions(::AbstractBGCFactory) = (;)

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

`explicit_override_keys` must list keys explicitly provided by the user (either
through `parameters` or `interaction_overrides`).

Resolution rules for each derived matrix key `K`:

1. If `K` is explicitly overridden, keep the provided value.
2. Else if `K` is missing, compute it.
3. Else if any declared dependency was explicitly overridden, recompute it.
4. Else keep the existing value.
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

"""Allometric palatability matrix derivation (`palatability_matrix[consumer, prey]`)."""
struct PalatabilityAllometric <: AbstractMatrixDeriver end

"""Binary assimilation-efficiency matrix derivation (`assimilation_matrix[consumer, prey]`)."""
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
