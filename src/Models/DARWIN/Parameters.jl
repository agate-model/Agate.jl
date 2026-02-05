"""Default parameter values for the DARWIN model.

This file defines `default_parameters(::DarwinFactory, ctx, FT)`, used by the
model-agnostic constructor.

Only keys required by the compiled DARWIN equations are used at runtime.

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
    resolve_diameter_indexed_vector,
    palatability_matrix_allometric_axes,
    assimilation_efficiency_matrix_binary_axes

import ...Utils: MatrixProvider, derived_matrix_specs

"""Parameter metadata for DARWIN.

The directory provides expected shapes and short descriptions for all parameters
required by the compiled DARWIN equations, plus a small set of interaction
traits used to derive the default interaction matrices.
"""
parameter_directory(::DarwinFactory) = (
    ParameterSpec(:DOC_remineralization, :scalar; doc="DOC remineralization rate."),
    ParameterSpec(:POC_remineralization, :scalar; doc="POC remineralization rate."),
    ParameterSpec(:DON_remineralization, :scalar; doc="DON remineralization rate."),
    ParameterSpec(:PON_remineralization, :scalar; doc="PON remineralization rate."),
    ParameterSpec(:DOP_remineralization, :scalar; doc="DOP remineralization rate."),
    ParameterSpec(:POP_remineralization, :scalar; doc="POP remineralization rate."),
    ParameterSpec(
        :nitrogen_to_carbon, :scalar; doc="Nitrogen-to-carbon stoichiometric ratio."
    ),
    ParameterSpec(
        :phosphorus_to_carbon, :scalar; doc="Phosphorus-to-carbon stoichiometric ratio."
    ),
    ParameterSpec(
        :DOM_POM_fractionation,
        :scalar;
        doc="Fraction of organic matter routed to DOM vs POM.",
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
        :half_saturation_DIN,
        :vector;
        doc="DIN half-saturation constant per plankton class.",
    ),
    ParameterSpec(
        :half_saturation_PO4,
        :vector;
        doc="PO4 half-saturation constant per plankton class.",
    ),
    ParameterSpec(
        :photosynthetic_slope,
        :vector;
        doc="Initial slope of the P-I curve per plankton class.",
    ),
    ParameterSpec(
        :chlorophyll_to_carbon_ratio,
        :vector;
        doc="Chlorophyll-to-carbon ratio per plankton class.",
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

# -----------------------------------------------------------------------------
# Derived interaction matrices
# -----------------------------------------------------------------------------

@inline function _derive_palatability_matrix(
    ::DarwinFactory, ctx::CommunityContext, params::NamedTuple
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
    ::DarwinFactory, ctx::CommunityContext, params::NamedTuple
)
    FT = ctx.FT
    return assimilation_efficiency_matrix_binary_axes(
        FT;
        assimilation_efficiency=params.assimilation_efficiency,
        consumer_indices=ctx.consumer_indices,
        prey_indices=ctx.prey_indices,
    )
end

function derived_matrix_specs(::DarwinFactory)
    return (;
        palatability_matrix=MatrixProvider(
            _derive_palatability_matrix;
            deps=(:optimum_predator_prey_ratio, :specificity, :protection),
        ),
        assimilation_matrix=MatrixProvider(
            _derive_assimilation_matrix; deps=(:assimilation_efficiency,)
        ),
    )
end

function default_parameters(::DarwinFactory, ctx::CommunityContext, ::Type{FT}) where {FT}
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
    linear_mortality = fill(FT(8e-7), ctx.n_total)
    quadratic_mortality = resolve_diameter_indexed_vector(
        FT, ctx.diameters, ctx.consumer_param_indices, FT(1e-6); default=FT(0)
    )

    # Phytoplankton growth (Geider light limitation + 2 nutrients)
    maximum_growth_rate = resolve_diameter_indexed_vector(
        FT,
        ctx.diameters,
        ctx.producer_param_indices,
        AllometricParam(PowerLaw(); prefactor=FT(2 / 86400), exponent=FT(-0.15));
        default=FT(0),
    )

    # Defaults to zero for non-phyto groups; MonodLimitation guards 0/0.
    half_saturation_DIN = resolve_diameter_indexed_vector(
        FT,
        ctx.diameters,
        ctx.producer_param_indices,
        AllometricParam(PowerLaw(); prefactor=FT(0.17), exponent=FT(0.27));
        default=FT(0),
    )
    half_saturation_PO4 = resolve_diameter_indexed_vector(
        FT,
        ctx.diameters,
        ctx.producer_param_indices,
        AllometricParam(PowerLaw(); prefactor=FT(0.17), exponent=FT(0.27));
        default=FT(0),
    )

    photosynthetic_slope = resolve_diameter_indexed_vector(
        FT, ctx.diameters, ctx.producer_param_indices, FT(0.1 / 86400); default=FT(0)
    )
    chlorophyll_to_carbon_ratio = resolve_diameter_indexed_vector(
        FT, ctx.diameters, ctx.producer_param_indices, FT(0.02); default=FT(0)
    )

    # Zooplankton grazing (preferential predation)
    maximum_predation_rate = resolve_diameter_indexed_vector(
        FT,
        ctx.diameters,
        ctx.consumer_param_indices,
        AllometricParam(PowerLaw(); prefactor=FT(30.84 / 86400), exponent=FT(-0.16));
        default=FT(0),
    )

    # Defaults to zero for non-zoo groups; HollingTypeII guards 0/0.
    holling_half_saturation = resolve_diameter_indexed_vector(
        FT,
        ctx.diameters,
        ctx.consumer_param_indices,
        AllometricParam(PowerLaw(); prefactor=FT(1.0), exponent=FT(-0.23));
        default=FT(0),
    )

    # ---------------------------------------------------------------------
    # Interaction matrices
    # ---------------------------------------------------------------------

    assimilation_efficiency = resolve_diameter_indexed_vector(
        FT, ctx.diameters, ctx.consumer_param_indices, FT(0.32); default=FT(0)
    )

    optimum_predator_prey_ratio = resolve_diameter_indexed_vector(
        FT, ctx.diameters, ctx.consumer_param_indices, FT(10.0); default=FT(0)
    )
    specificity = resolve_diameter_indexed_vector(
        FT, ctx.diameters, ctx.consumer_param_indices, FT(0.3); default=FT(0)
    )
    protection = fill(FT(0), ctx.n_total)

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
        optimum_predator_prey_ratio,
        specificity,
        protection,
        assimilation_efficiency,
    )
end
