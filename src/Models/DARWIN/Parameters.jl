"""Default parameter values for the DARWIN model.

This file defines `default_parameters(::DarwinFactory, ctx, FT)`, used by the
model-agnostic constructor.

Only keys required by the compiled DARWIN equations are used at runtime.

In addition, this file exposes a small set of *interaction traits* (vectors)
that may be overridden to regenerate the interaction matrices during
construction.
"""

import ...Constructor: default_parameters
import ...FactoryInterface: parameter_directory, ParameterSpec
using ...Utils: InteractionContext
using ...Library.Allometry:
    AllometricParam,
    PowerLaw,
    resolve_param,
    palatability_matrix_allometric_axes,
    assimilation_efficiency_matrix_binary_axes

import ...Utils: MatrixFn, derived_matrix_specs

"""Parameter metadata for DARWIN.

The directory provides expected shapes and short descriptions for all parameters
required by the compiled DARWIN equations, plus a small set of interaction
traits used to derive the default interaction matrices.
"""
parameter_directory(::DarwinFactory) = (
    ParameterSpec(
        :DOC_remineralization, :scalar; kind=:real, doc="DOC remineralization rate."
    ),
    ParameterSpec(
        :POC_remineralization, :scalar; kind=:real, doc="POC remineralization rate."
    ),
    ParameterSpec(
        :DON_remineralization, :scalar; kind=:real, doc="DON remineralization rate."
    ),
    ParameterSpec(
        :PON_remineralization, :scalar; kind=:real, doc="PON remineralization rate."
    ),
    ParameterSpec(
        :DOP_remineralization, :scalar; kind=:real, doc="DOP remineralization rate."
    ),
    ParameterSpec(
        :POP_remineralization, :scalar; kind=:real, doc="POP remineralization rate."
    ),
    ParameterSpec(
        :nitrogen_to_carbon,
        :scalar;
        kind=:real,
        doc="Nitrogen-to-carbon stoichiometric ratio.",
    ),
    ParameterSpec(
        :phosphorus_to_carbon,
        :scalar;
        kind=:real,
        doc="Phosphorus-to-carbon stoichiometric ratio.",
    ),
    ParameterSpec(
        :DOM_POM_fractionation,
        :scalar;
        kind=:real,
        doc="Fraction of organic matter routed to DOM vs POM.",
    ),
    ParameterSpec(
        :linear_mortality,
        :vector;
        kind=:real,
        doc="Linear mortality coefficient per plankton class.",
    ),
    ParameterSpec(
        :quadratic_mortality,
        :vector;
        kind=:real,
        doc="Quadratic mortality coefficient per plankton class.",
    ),
    ParameterSpec(
        :maximum_growth_rate,
        :vector;
        kind=:real,
        doc="Maximum phytoplankton growth rate per plankton class.",
    ),
    ParameterSpec(
        :half_saturation_DIN,
        :vector;
        kind=:real,
        doc="DIN half-saturation constant per plankton class.",
    ),
    ParameterSpec(
        :half_saturation_PO4,
        :vector;
        kind=:real,
        doc="PO4 half-saturation constant per plankton class.",
    ),
    ParameterSpec(
        :photosynthetic_slope,
        :vector;
        kind=:real,
        doc="Initial slope of the P-I curve per plankton class.",
    ),
    ParameterSpec(
        :chlorophyll_to_carbon_ratio,
        :vector;
        kind=:real,
        doc="Chlorophyll-to-carbon ratio per plankton class.",
    ),
    ParameterSpec(
        :maximum_predation_rate,
        :vector;
        kind=:real,
        doc="Maximum zooplankton grazing rate per plankton class.",
    ),
    ParameterSpec(
        :holling_half_saturation,
        :vector;
        kind=:real,
        doc="Holling type II half-saturation constant per plankton class.",
    ),
    ParameterSpec(
        :palatability_matrix,
        :matrix;
        kind=:real,
        axes=(:consumer, :prey),
        doc="Preference of each consumer for each prey class.",
    ),
    ParameterSpec(
        :assimilation_matrix,
        :matrix;
        kind=:real,
        axes=(:consumer, :prey),
        doc="Assimilation efficiency of each consumer on each prey class.",
    ),

    # Interaction traits used to derive the default matrices.
    ParameterSpec(
        :can_eat,
        :vector;
        kind=:bool,
        doc="Whether each plankton class can act as a consumer in grazing interactions.",
    ),
    ParameterSpec(
        :can_be_eaten,
        :vector;
        kind=:bool,
        doc="Whether each plankton class can be grazed (acts as prey).",
    ),
    ParameterSpec(
        :optimum_predator_prey_ratio,
        :vector;
        kind=:real,
        doc="Preferred predator:prey diameter ratio per consumer (used to derive palatability_matrix).",
    ),
    ParameterSpec(
        :specificity,
        :vector;
        kind=:real,
        doc="Unimodal palatability specificity per consumer (used to derive palatability_matrix).",
    ),
    ParameterSpec(
        :protection,
        :vector;
        kind=:real,
        doc="Prey protection factor (used to derive palatability_matrix).",
    ),
    ParameterSpec(
        :assimilation_efficiency,
        :vector;
        kind=:real,
        doc="Assimilation efficiency per consumer (used to derive assimilation_matrix).",
    ),
)

@inline function _group_value(group_map::NamedTuple, group::Symbol, default)
    return Base.hasproperty(group_map, group) ? getproperty(group_map, group) : default
end

# -----------------------------------------------------------------------------
# Derived interaction matrices
# -----------------------------------------------------------------------------

@inline function _derive_palatability_matrix(
    ::DarwinFactory, ctx::InteractionContext, params::NamedTuple
)
    FT = ctx.FT
    return palatability_matrix_allometric_axes(
        FT,
        ctx.diameters;
        can_eat=params.can_eat,
        can_be_eaten=params.can_be_eaten,
        optimum_predator_prey_ratio=params.optimum_predator_prey_ratio,
        specificity=params.specificity,
        protection=params.protection,
        consumer_indices=ctx.consumer_indices,
        prey_indices=ctx.prey_indices,
    )
end

@inline function _derive_assimilation_matrix(
    ::DarwinFactory, ctx::InteractionContext, params::NamedTuple
)
    FT = ctx.FT
    return assimilation_efficiency_matrix_binary_axes(
        FT;
        can_eat=params.can_eat,
        can_be_eaten=params.can_be_eaten,
        assimilation_efficiency=params.assimilation_efficiency,
        consumer_indices=ctx.consumer_indices,
        prey_indices=ctx.prey_indices,
    )
end

function derived_matrix_specs(::DarwinFactory)
    return (;
        palatability_matrix=MatrixFn(
            _derive_palatability_matrix;
            deps=(
                :can_eat,
                :can_be_eaten,
                :optimum_predator_prey_ratio,
                :specificity,
                :protection,
            ),
        ),
        assimilation_matrix=MatrixFn(
            _derive_assimilation_matrix;
            deps=(:can_eat, :can_be_eaten, :assimilation_efficiency),
        ),
    )
end

function _resolve_groupvec(
    ::Type{FT}, ctx::InteractionContext, group_map::NamedTuple; default
) where {FT}
    n = ctx.n_total
    out = Vector{FT}(undef, n)
    @inbounds for i in 1:n
        g = ctx.group_symbols[i]
        v = _group_value(group_map, g, default)
        out[i] = resolve_param(FT, v, ctx.diameters[i])
    end
    return out
end

function _resolve_groupvec_bool(
    ctx::InteractionContext, group_map::NamedTuple; default=false
)
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

    optimum_predator_prey_ratio = _resolve_groupvec(FT, ctx, (; Z=FT(10.0)); default=0.0)
    specificity = _resolve_groupvec(FT, ctx, (; Z=FT(0.3)); default=0.0)
    protection = _resolve_groupvec(FT, ctx, (;); default=0.0)

    # Canonical interaction storage is consumer-by-prey.
    consumer_idx = ctx.consumer_indices
    prey_idx = ctx.prey_indices

    palatability_matrix = palatability_matrix_allometric_axes(
        FT,
        ctx.diameters;
        can_eat=can_eat,
        can_be_eaten=can_be_eaten,
        optimum_predator_prey_ratio=optimum_predator_prey_ratio,
        specificity=specificity,
        protection=protection,
        consumer_indices=consumer_idx,
        prey_indices=prey_idx,
    )

    assimilation_matrix = assimilation_efficiency_matrix_binary_axes(
        FT;
        can_eat=can_eat,
        can_be_eaten=can_be_eaten,
        assimilation_efficiency=assimilation_efficiency,
        consumer_indices=consumer_idx,
        prey_indices=prey_idx,
    )

    return (;
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
        can_eat,
        can_be_eaten,
        optimum_predator_prey_ratio,
        specificity,
        protection,
        assimilation_efficiency,
    )
end
