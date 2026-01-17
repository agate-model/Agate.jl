"""Shared default interaction-matrix providers.

These providers are used as *defaults* for model parameter registries.
"""
module InteractionDefaults

using ...Parameters: MatrixFn
using ...Library.Allometry:
    allometric_palatability_unimodal_protection,
    palatability_matrix_allometric,
    assimilation_efficiency_matrix_binary

@inline function _default_palatability_matrix(ctx, depvals)
    can_eat, optimum_predator_prey_ratio, specificity, protection = depvals
    return palatability_matrix_allometric(ctx.FT, ctx.diameters;
        can_eat=can_eat,
        optimum_predator_prey_ratio=optimum_predator_prey_ratio,
        specificity=specificity,
        protection=protection,
        palatability_fn=allometric_palatability_unimodal_protection,
    )
end

@inline function _default_assimilation_matrix(ctx, depvals)
    can_eat, can_be_eaten, assimilation_efficiency = depvals
    return assimilation_efficiency_matrix_binary(eltype(ctx.diameters);
        can_eat=can_eat,
        can_be_eaten=can_be_eaten,
        assimilation_efficiency=assimilation_efficiency,
    )
end

const DEFAULT_PALATABILITY_PROVIDER = MatrixFn(
    _default_palatability_matrix;
    deps=(:can_eat, :optimum_predator_prey_ratio, :specificity, :protection),
)

const DEFAULT_ASSIMILATION_PROVIDER = MatrixFn(
    _default_assimilation_matrix;
    deps=(:can_eat, :can_be_eaten, :assimilation_efficiency),
)

"""Return the canonical default palatability matrix provider."""
default_palatability_provider() = DEFAULT_PALATABILITY_PROVIDER

"""Return the canonical default assimilation-efficiency matrix provider."""
default_assimilation_provider() = DEFAULT_ASSIMILATION_PROVIDER

export default_palatability_provider, default_assimilation_provider

end # module InteractionDefaults
