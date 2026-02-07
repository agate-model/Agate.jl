"""Shared helpers to derive role-aware interaction matrices from trait vectors.

These helpers implement the common pattern used by several plankton ecosystem
models in Agate:

- expose small trait vectors (e.g. predator:prey size ratio preferences)
- derive full consumer-by-prey interaction matrices during construction

The returned matrices are host arrays. They can later be moved to the target
architecture via `Adapt` during model construction.
"""

import ..Utils: AbstractBGCFactory, CommunityContext

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
