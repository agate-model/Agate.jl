"""Default parameter values for the NiPiZD model.

The public constructor (`NiPiZD.construct`) takes a `parameters::NamedTuple` override.
This file provides a single method that generates a complete, resolved parameter
`NamedTuple` for a particular community (sizes + groups) and floating-point type.

Only parameters required by the model's compiled equations are produced. Interaction
traits used to derive default matrices are kept internal to this method.
"""

import ...Constructor: default_parameters
using ...Utils: InteractionContext
using ...Library.Allometry:
    AllometricParam,
    PowerLaw,
    resolve_param,
    palatability_matrix_allometric,
    assimilation_efficiency_matrix_binary

@inline function _group_value(group_map::NamedTuple, group::Symbol, default)
    return Base.hasproperty(group_map, group) ? getproperty(group_map, group) : default
end

function _resolve_groupvec(::Type{FT}, ctx::InteractionContext, group_map::NamedTuple; default=0.0) where {FT}
    n = ctx.n_total
    out = Vector{FT}(undef, n)
    for i in 1:n
        g = ctx.group_symbols[i]
        v = _group_value(group_map, g, default)
        out[i] = resolve_param(FT, v, ctx.diameters[i])
    end
    return out
end

function _resolve_groupvec_bool(ctx::InteractionContext, group_map::NamedTuple; default=false)
    n = ctx.n_total
    out = Vector{Bool}(undef, n)
    for i in 1:n
        g = ctx.group_symbols[i]
        out[i] = Bool(_group_value(group_map, g, default))
    end
    return out
end

function default_parameters(::NiPiZDFactory, ctx::InteractionContext, ::Type{FT}) where {FT}
    groups = (:Z, :P)

    detritus_remineralization = FT(0.1213 / 86400)
    mortality_export_fraction = FT(0.2)

    linear_mortality = _resolve_groupvec(FT, ctx, (; Z=8e-7, P=8e-7); default=0.0)
    quadratic_mortality = _resolve_groupvec(FT, ctx, (; Z=1e-6, P=0.0); default=0.0)

    maximum_growth_rate = _resolve_groupvec(
        FT,
        ctx,
        (; P=AllometricParam(PowerLaw(); prefactor=2 / 86400, exponent=-0.15), Z=0.0);
        default=0.0,
    )
    nutrient_half_saturation = _resolve_groupvec(
        FT,
        ctx,
        (; P=AllometricParam(PowerLaw(); prefactor=0.17, exponent=0.27), Z=0.0);
        default=0.0,
    )
    alpha = _resolve_groupvec(FT, ctx, (; P=0.1953 / 86400, Z=0.0); default=0.0)

    maximum_predation_rate = _resolve_groupvec(
        FT,
        ctx,
        (; Z=AllometricParam(PowerLaw(); prefactor=30.84 / 86400, exponent=-0.16), P=0.0);
        default=0.0,
    )
    holling_half_saturation = _resolve_groupvec(FT, ctx, (; Z=5.0, P=0.0); default=0.0)

    # --- Default interaction matrices (internal trait defaults) ------------
    can_eat = _resolve_groupvec_bool(ctx, (; Z=true, P=false); default=false)
    can_be_eaten = _resolve_groupvec_bool(ctx, (; Z=false, P=true); default=false)

    optimum_ratio = _resolve_groupvec(FT, ctx, (; Z=10.0, P=0.0); default=0.0)
    specificity = _resolve_groupvec(FT, ctx, (; Z=0.3, P=0.0); default=0.0)
    protection = _resolve_groupvec(FT, ctx, (; Z=1.0, P=0.0); default=0.0)
    assimilation_efficiency = _resolve_groupvec(FT, ctx, (; Z=0.32, P=0.0); default=0.0)

    palatability_matrix = palatability_matrix_allometric(
        FT,
        ctx.diameters;
        can_eat=can_eat,
        optimum_predator_prey_ratio=optimum_ratio,
        specificity=specificity,
        protection=protection,
    )

    assimilation_matrix = assimilation_efficiency_matrix_binary(
        can_eat,
        can_be_eaten,
        assimilation_efficiency,
    )

    return (
        ;
        detritus_remineralization,
        mortality_export_fraction,
        linear_mortality,
        quadratic_mortality,
        maximum_growth_rate,
        nutrient_half_saturation,
        alpha,
        maximum_predation_rate,
        holling_half_saturation,
        palatability_matrix,
        assimilation_matrix,
    )
end
