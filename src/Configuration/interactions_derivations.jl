# """Utilities for *derived* interaction matrices.
# 
# Some models expose low-level interaction *traits* (e.g. predator-prey size ratio
# preferences) and derive interaction matrices (palatability, assimilation) from
# those traits.
# 
# Derived matrices are computed during construction when missing. They are
# recomputed when:
# 
# - a matrix is *not* explicitly overridden, and
# - at least one of its declared dependencies *is* explicitly overridden.
# 
# This allows users to tweak traits (small surface, easy to reason about) without
# needing to hand-build full matrices.
# """

# """Return the dependency keys for a derivation callable.
# 
# Derivation dependencies are used to decide when a derived matrix must be
# recomputed if its inputs are explicitly overridden.
# 
# Override this for specific derivation functions.
# """
derivation_deps(::Function) = ()

"""Return a `NamedTuple` mapping matrix keys to derivation callables.

Factories override this to declare which matrices can be derived from other
parameters.
"""
derived_matrix_specs(::AbstractBGCFactory) = (;)

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
    derived_matrix_providers = derived_matrix_specs(factory)
    isempty(derived_matrix_providers) && return params

    override_set = Set(explicit_override_keys)
    matrix_keys = propertynames(derived_matrix_providers)

    vals = ntuple(
        i -> begin
            key = matrix_keys[i]
            provider = getproperty(derived_matrix_providers, key)
            deps = derivation_deps(provider)
            deps_norm = if deps === nothing
                ()
            elseif deps isa Symbol
                (deps,)
            elseif deps isa AbstractVector
                Tuple(deps)
            else
                deps
            end

            # Rule 1: explicit override wins.
            if key in override_set && hasproperty(params, key)
                return getproperty(params, key)
            end

            # Rule 2: compute missing matrices.
            if !hasproperty(params, key)
                return provider(factory, context, params)
            end

            # Rule 3: recompute when any dependency was explicitly overridden.
            if !isempty(deps_norm) && any(d -> d in override_set, deps_norm)
                return provider(factory, context, params)
            end

            # Rule 4: otherwise keep the current value.
            return getproperty(params, key)
        end,
        length(matrix_keys),
    )

    updates = NamedTuple{matrix_keys}(vals)
    return merge(params, updates)
end

# """Shared helpers to derive role-aware interaction matrices from trait vectors.
# 
# These helpers implement the common pattern used by several plankton ecosystem
# models in Agate:
# 
# - expose small trait vectors (e.g. predator:prey size ratio preferences)
# - derive full consumer-by-prey interaction matrices during construction
# 
# The returned matrices are host arrays. They can later be moved to the target
# architecture via `Adapt` during model construction.
# """

using ..Library.Allometry:
    palatability_matrix_allometric_axes, assimilation_efficiency_matrix_binary_axes

"""Derive `palatability_matrix` from trait vectors in `params`.

Required fields in `params`:
- `optimum_predator_prey_ratio`
- `specificity`
- `protection`
"""
@inline function derive_palatability_matrix_allometric(
    ::AbstractBGCFactory, community_context::CommunityContext, params::NamedTuple
)
    FT = community_context.FT
    return palatability_matrix_allometric_axes(
        FT,
        community_context.diameters;
        optimum_predator_prey_ratio=params.optimum_predator_prey_ratio,
        specificity=params.specificity,
        protection=params.protection,
        consumer_indices=community_context.consumer_indices,
        prey_indices=community_context.prey_indices,
    )
end

derivation_deps(::typeof(derive_palatability_matrix_allometric)) =
    (:optimum_predator_prey_ratio, :specificity, :protection)

"""Derive `assimilation_matrix` from trait vectors in `params`.

Required fields in `params`:
- `assimilation_efficiency`
"""
@inline function derive_assimilation_matrix_binary(
    ::AbstractBGCFactory, community_context::CommunityContext, params::NamedTuple
)
    FT = community_context.FT
    return assimilation_efficiency_matrix_binary_axes(
        FT;
        assimilation_efficiency=params.assimilation_efficiency,
        consumer_indices=community_context.consumer_indices,
        prey_indices=community_context.prey_indices,
    )
end

derivation_deps(::typeof(derive_assimilation_matrix_binary)) = (:assimilation_efficiency,)
