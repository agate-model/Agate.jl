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
import ...Interface: parameter_directory, ParameterSpec
using ...Utils: CommunityContext
using ...Library.Allometry:
    AllometricParam,
    PowerLaw,
    resolve_param,
    palatability_matrix_allometric_axes,
    assimilation_efficiency_matrix_binary_axes

import ...Utils: MatrixProvider, derived_matrix_specs

function _resolve_allvec(::Type{FT}, ctx::CommunityContext, value) where {FT}
    n = ctx.n_total
    out = Vector{FT}(undef, n)
    @inbounds for i in 1:n
        out[i] = resolve_param(FT, value, ctx.diameters[i])
    end
    return out
end

function _resolve_indexedvec(
    ::Type{FT}, ctx::CommunityContext, indices::AbstractVector{Int}, value; default::FT
) where {FT}
    n = ctx.n_total
    out = Vector{FT}(undef, n)
    @inbounds for i in 1:n
        out[i] = default
    end
    @inbounds for i in indices
        out[i] = resolve_param(FT, value, ctx.diameters[i])
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
        doc="Detritus remineralization rate.",
    ),
    ParameterSpec(
        :mortality_export_fraction,
        :scalar;
        doc="Fraction of mortality routed to detritus export.",
    ),
    ParameterSpec(
        :linear_mortality,
        :vector;
        doc="Linear mortality coefficient per plankton class.",
    ),
    ParameterSpec(
        :quadratic_mortality,
        :vector;
        doc="Quadratic mortality coefficient per plankton class.",
    ),
    ParameterSpec(
        :maximum_growth_rate,
        :vector;
        doc="Maximum phytoplankton growth rate per plankton class.",
    ),
    ParameterSpec(
        :nutrient_half_saturation,
        :vector;
        doc="Nutrient half-saturation constant per plankton class.",
    ),
    ParameterSpec(
        :alpha,
        :vector;
        doc="Initial slope of the P-I curve per plankton class.",
    ),
    ParameterSpec(
        :maximum_predation_rate,
        :vector;
        doc="Maximum zooplankton grazing rate per plankton class.",
    ),
    ParameterSpec(
        :holling_half_saturation,
        :vector;
        doc="Holling type II half-saturation constant per plankton class.",
    ),
    ParameterSpec(
        :palatability_matrix,
        :matrix;
        axes=(:consumer, :prey),
        doc="Preference of each consumer for each prey class.",
    ),
    ParameterSpec(
        :assimilation_matrix,
        :matrix;
        axes=(:consumer, :prey),
        doc="Assimilation efficiency of each consumer on each prey class.",
    ),

    # Traits used to derive the default interaction matrices.
    ParameterSpec(
        :optimum_predator_prey_ratio,
        :vector;
        doc="Preferred predator:prey diameter ratio per consumer (used to derive palatability_matrix).",
    ),
    ParameterSpec(
        :specificity,
        :vector;
        doc="Unimodal palatability specificity per consumer (used to derive palatability_matrix).",
    ),
    ParameterSpec(
        :protection,
        :vector;
        doc="Prey protection factor (used to derive palatability_matrix).",
    ),
    ParameterSpec(
        :assimilation_efficiency,
        :vector;
        doc="Assimilation efficiency per consumer (used to derive assimilation_matrix).",
    ),
)

function default_parameters(::NiPiZDFactory, ctx::CommunityContext, ::Type{FT}) where {FT}
    detritus_remineralization = FT(0.1213 / 86400)
    mortality_export_fraction = FT(0.2)

    # Group vectors (length = ctx.n_total)
    linear_mortality = fill(FT(8e-7), ctx.n_total)
    quadratic_mortality = _resolve_indexedvec(FT, ctx, ctx.consumer_param_indices, FT(1e-6); default=FT(0))

    maximum_growth_rate = _resolve_indexedvec(
        FT,
        ctx,
        ctx.producer_param_indices,
        AllometricParam(PowerLaw(); prefactor=FT(2 / 86400), exponent=FT(-0.15));
        default=FT(0),
    )

    nutrient_half_saturation = _resolve_indexedvec(
        FT,
        ctx,
        ctx.producer_param_indices,
        AllometricParam(PowerLaw(); prefactor=FT(0.17), exponent=FT(0.27));
        default=FT(0),
    )

    alpha = _resolve_indexedvec(FT, ctx, ctx.producer_param_indices, FT(0.1953 / 86400); default=FT(0))

    maximum_predation_rate = _resolve_indexedvec(
        FT,
        ctx,
        ctx.consumer_param_indices,
        AllometricParam(PowerLaw(); prefactor=FT(30.84 / 86400), exponent=FT(-0.16));
        default=FT(0),
    )

    holling_half_saturation = _resolve_indexedvec(FT, ctx, ctx.consumer_param_indices, FT(5.0); default=FT(0))

    # --- Default interaction matrices (internal trait defaults) ------------

    optimum_predator_prey_ratio = _resolve_indexedvec(
        FT, ctx, ctx.consumer_param_indices, FT(10.0); default=FT(0)
    )
    specificity = _resolve_indexedvec(FT, ctx, ctx.consumer_param_indices, FT(0.3); default=FT(0))
    protection = _resolve_indexedvec(FT, ctx, ctx.consumer_param_indices, FT(1.0); default=FT(0))
    assimilation_efficiency = _resolve_indexedvec(
        FT, ctx, ctx.consumer_param_indices, FT(0.32); default=FT(0)
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
    ::NiPiZDFactory, ctx::CommunityContext, params::NamedTuple
)
    FT = ctx.FT
    return palatability_matrix_allometric_axes(
        FT,
        ctx.diameters;

        optimum_predator_prey_ratio=params.optimum_predator_prey_ratio,
        specificity=params.specificity,
        protection=params.protection,
        consumer_indices=ctx.consumer_indices,
        prey_indices=ctx.prey_indices,
    )
end

@inline function _derive_assimilation_matrix(
    ::NiPiZDFactory, ctx::CommunityContext, params::NamedTuple
)
    FT = ctx.FT
    return assimilation_efficiency_matrix_binary_axes(
        FT;

        assimilation_efficiency=params.assimilation_efficiency,
        consumer_indices=ctx.consumer_indices,
        prey_indices=ctx.prey_indices,
    )
end

function derived_matrix_specs(::NiPiZDFactory)
    return (;
        palatability_matrix=MatrixProvider(
            _derive_palatability_matrix;
            deps=(:optimum_predator_prey_ratio, :specificity, :protection),
        ),
        assimilation_matrix=MatrixProvider(
            _derive_assimilation_matrix;
            deps=(:assimilation_efficiency,),
        ),
    )
end
