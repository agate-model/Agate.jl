"""Default parameter values for the DARWIN model.

This file defines `default_parameters(::DarwinFactory, ctx, FT)`, used by the
model-agnostic constructor.

Only keys required by the compiled DARWIN equations are provided.
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

function _resolve_groupvec(::Type{FT}, ctx::InteractionContext, group_map::NamedTuple; default) where {FT}
    n = ctx.n_total
    out = Vector{FT}(undef, n)
    @inbounds for i in 1:n
        g = ctx.group_symbols[i]
        v = _group_value(group_map, g, default)
        out[i] = resolve_param(FT, v, ctx.diameters[i])
    end
    return out
end

function _resolve_groupvec_bool(ctx::InteractionContext, group_map::NamedTuple; default=false)
    n = ctx.n_total
    out = Vector{Bool}(undef, n)
    @inbounds for i in 1:n
        g = ctx.group_symbols[i]
        out[i] = Bool(_group_value(group_map, g, default))
    end
    return out
end

function default_parameters(::DarwinFactory, ctx::InteractionContext, ::Type{FT}) where {FT}
    # ---------------------------------------------------------------------
    # Scalars
    # ---------------------------------------------------------------------

    detritus_remin = FT(0.1213 / 86400)

    DOC_remineralization = detritus_remin
    POC_remineralization = detritus_remin
    DON_remineralization = detritus_remin
    PON_remineralization = detritus_remin
    DOP_remineralization = detritus_remin
    POP_remineralization = detritus_remin

    nitrogen_to_carbon = FT(16 / 106)
    phosphorus_to_carbon = FT(1 / 106)

    DOM_POM_fractionation = FT(0.5)

    # ---------------------------------------------------------------------
    # Group vectors (length = ctx.n_total)
    # ---------------------------------------------------------------------

    # Mortality
    linear_mortality = _resolve_groupvec(FT, ctx, (; Z=8e-7, P=8e-7); default=0.0)
    quadratic_mortality = _resolve_groupvec(FT, ctx, (; Z=1e-6); default=0.0)

    # Phytoplankton growth (Geider light limitation + 2 nutrients)
    maximum_growth_rate = _resolve_groupvec(
        FT,
        ctx,
        (; P=AllometricParam(PowerLaw(); prefactor=FT(2 / 86400), exponent=FT(-0.15)));
        default=0.0,
    )

    # Defaults to zero for non-phyto groups; MonodLimitation guards 0/0.
    half_saturation_DIN = _resolve_groupvec(
        FT,
        ctx,
        (; P=AllometricParam(PowerLaw(); prefactor=FT(0.17), exponent=FT(0.27)));
        default=FT(0.0),
    )
    half_saturation_PO4 = _resolve_groupvec(
        FT,
        ctx,
        (; P=AllometricParam(PowerLaw(); prefactor=FT(0.17), exponent=FT(0.27)));
        default=FT(0.0),
    )

    photosynthetic_slope = _resolve_groupvec(FT, ctx, (; P=FT(0.1 / 86400)); default=0.0)
    chlorophyll_to_carbon_ratio = _resolve_groupvec(FT, ctx, (; P=FT(0.02)); default=0.0)

    # Zooplankton grazing (preferential predation)
    maximum_predation_rate = _resolve_groupvec(
        FT,
        ctx,
        (; Z=AllometricParam(PowerLaw(); prefactor=FT(30.84 / 86400), exponent=FT(-0.16)));
        default=0.0,
    )

    # Defaults to zero for non-zoo groups; HollingTypeII guards 0/0.
    holling_half_saturation = _resolve_groupvec(
        FT,
        ctx,
        (; Z=AllometricParam(PowerLaw(); prefactor=FT(1.0), exponent=FT(-0.23)));
        default=FT(0.0),
    )

    # ---------------------------------------------------------------------
    # Interaction matrices
    # ---------------------------------------------------------------------

    can_eat = _resolve_groupvec_bool(ctx, (; Z=true); default=false)
    can_be_eaten = _resolve_groupvec_bool(ctx, (; P=true); default=false)

    assimilation_efficiency = _resolve_groupvec(FT, ctx, (; Z=FT(0.32)); default=0.0)

    optimum_ratio = _resolve_groupvec(FT, ctx, (; Z=FT(10.0)); default=0.0)
    specificity = _resolve_groupvec(FT, ctx, (; Z=FT(0.3)); default=0.0)
    protection = _resolve_groupvec(FT, ctx, (;); default=0.0)

    palatability_matrix = palatability_matrix_allometric(
        FT,
        ctx.diameters;
        can_eat=can_eat,
        optimum_predator_prey_ratio=optimum_ratio,
        specificity=specificity,
        protection=protection,
    )

    assimilation_matrix = assimilation_efficiency_matrix_binary(can_eat, can_be_eaten, assimilation_efficiency)

    return (
        ;
        DOC_remineralization,
        POC_remineralization,
        DON_remineralization,
        PON_remineralization,
        DOP_remineralization,
        POP_remineralization,
        nitrogen_to_carbon,
        phosphorus_to_carbon,
        DOM_POM_fractionation,
        linear_mortality,
        quadratic_mortality,
        maximum_growth_rate,
        half_saturation_DIN,
        half_saturation_PO4,
        photosynthetic_slope,
        chlorophyll_to_carbon_ratio,
        maximum_predation_rate,
        holling_half_saturation,
        palatability_matrix,
        assimilation_matrix,
    )
end
