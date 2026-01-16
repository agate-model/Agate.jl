# This file is included into the `Agate.Parameters` module.

# ----------------------------------------------------------------------------
# Default derived interaction matrices
# ----------------------------------------------------------------------------

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
    FT = eltype(ctx.diameters)
    return assimilation_efficiency_matrix_binary(FT;
        can_eat=can_eat,
        can_be_eaten=can_be_eaten,
        assimilation_efficiency=assimilation_efficiency,
    )
end

"""Return the canonical default palatability matrix provider."""
default_palatability_provider() =
    MatrixFn(_default_palatability_matrix; deps=[:can_eat, :optimum_predator_prey_ratio, :specificity, :protection])

"""Return the canonical default assimilation efficiency matrix provider."""
default_assimilation_provider() =
    MatrixFn(_default_assimilation_matrix; deps=[:can_eat, :can_be_eaten, :assimilation_efficiency])
