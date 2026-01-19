"""Default parameter values for the NiPiZD model.

The public constructor (`NiPiZD.construct`) takes a `parameters::NamedTuple` override.
This file provides a single method that generates a complete, resolved parameter
`NamedTuple` for a particular community (sizes + groups) and floating-point type.

Only parameters required by the model's compiled equations are used at runtime.

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

@inline function _group_value(group_map::NamedTuple, group::Symbol, default)
    return Base.hasproperty(group_map, group) ? getproperty(group_map, group) : default
end

function _resolve_groupvec(
    ::Type{FT}, ctx::InteractionContext, group_map::NamedTuple; default=0.0
) where {FT}
    n = ctx.n_total
    out = Vector{FT}(undef, n)
    for i in 1:n
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
    for i in 1:n
        g = ctx.group_symbols[i]
        out[i] = Bool(_group_value(group_map, g, default))
    end
    return out
end

"""Parameter metadata for NiPiZD.

The directory provides expected shapes and short descriptions for all parameters
required by the compiled NiPiZD equations, plus a small set of interaction
traits used to derive the default interaction matrices.
"""
parameter_directory(::NiPiZDFactory) = (
    ParameterSpec(
        :detritus_remineralization,
        :scalar;
        kind=:real,
        doc="Detritus remineralization rate.",
    ),
    ParameterSpec(
        :mortality_export_fraction,
        :scalar;
        kind=:real,
        doc="Fraction of mortality routed to detritus export.",
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
        :nutrient_half_saturation,
        :vector;
        kind=:real,
        doc="Nutrient half-saturation constant per plankton class.",
    ),
    ParameterSpec(
        :alpha,
        :vector;
        kind=:real,
        doc="Initial slope of the P-I curve per plankton class.",
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

function default_parameters(::NiPiZDFactory, ctx::InteractionContext, ::Type{FT}) where {FT}
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

    optimum_predator_prey_ratio = _resolve_groupvec(FT, ctx, (; Z=10.0, P=0.0); default=0.0)
    specificity = _resolve_groupvec(FT, ctx, (; Z=0.3, P=0.0); default=0.0)
    protection = _resolve_groupvec(FT, ctx, (; Z=1.0, P=0.0); default=0.0)
    assimilation_efficiency = _resolve_groupvec(FT, ctx, (; Z=0.32, P=0.0); default=0.0)

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
        can_eat,
        can_be_eaten,
        optimum_predator_prey_ratio,
        specificity,
        protection,
        assimilation_efficiency,
    )
end

# -----------------------------------------------------------------------------
# Derived interaction matrices
# -----------------------------------------------------------------------------

@inline function _derive_palatability_matrix(
    ::NiPiZDFactory, ctx::InteractionContext, params::NamedTuple
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
    ::NiPiZDFactory, ctx::InteractionContext, params::NamedTuple
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

function derived_matrix_specs(::NiPiZDFactory)
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
