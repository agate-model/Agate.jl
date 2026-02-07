"""Default parameter values for the NiPiZD model.

The public constructor (`NiPiZD.construct`) accepts a `parameters::NamedTuple` override.
This file registers constructor-time default parameters via `parameter_default_registry`.
Defaults are computed for a particular community (sizes + groups) and floating-point type.

Only parameters required by the model's compiled equations are used at runtime.

In addition, this file exposes a small set of *interaction traits* (vectors)
that may be overridden to regenerate the interaction matrices during
construction.
"""

import ...Constructor: parameter_default_registry, FnDefault, NoDefault
import ...Interface: parameter_directory, ParameterSpec
using ...Utils: CommunityContext
using ...Library.Allometry:
    AllometricParam,
    PowerLaw,
    resolve_diameter_indexed_vector,
    palatability_matrix_allometric_axes,
    assimilation_efficiency_matrix_binary_axes

import ...Utils: MatrixProvider, derived_matrix_specs

"""Parameter metadata for NiPiZD.

The directory provides expected shapes and short descriptions for all parameters
required by the compiled NiPiZD equations, plus a small set of interaction
traits used to derive the default interaction matrices.
"""
parameter_directory(::NiPiZDFactory) = (
    ParameterSpec(
        :detritus_remineralization, :scalar; doc="Detritus remineralization rate."
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
        :alpha, :vector; doc="Initial slope of the P-I curve per plankton class."
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

function parameter_default_registry(factory::NiPiZDFactory)
    return (
        detritus_remineralization = FnDefault((factory, community_context, FT, cache) ->
            nipizd_parameter_default_values(factory, community_context, FT, cache).detritus_remineralization
        ),
        mortality_export_fraction = FnDefault((factory, community_context, FT, cache) ->
            nipizd_parameter_default_values(factory, community_context, FT, cache).mortality_export_fraction
        ),
        linear_mortality = FnDefault((factory, community_context, FT, cache) ->
            nipizd_parameter_default_values(factory, community_context, FT, cache).linear_mortality
        ),
        quadratic_mortality = FnDefault((factory, community_context, FT, cache) ->
            nipizd_parameter_default_values(factory, community_context, FT, cache).quadratic_mortality
        ),
        maximum_growth_rate = FnDefault((factory, community_context, FT, cache) ->
            nipizd_parameter_default_values(factory, community_context, FT, cache).maximum_growth_rate
        ),
        nutrient_half_saturation = FnDefault((factory, community_context, FT, cache) ->
            nipizd_parameter_default_values(factory, community_context, FT, cache).nutrient_half_saturation
        ),
        alpha = FnDefault((factory, community_context, FT, cache) ->
            nipizd_parameter_default_values(factory, community_context, FT, cache).alpha
        ),
        maximum_predation_rate = FnDefault((factory, community_context, FT, cache) ->
            nipizd_parameter_default_values(factory, community_context, FT, cache).maximum_predation_rate
        ),
        holling_half_saturation = FnDefault((factory, community_context, FT, cache) ->
            nipizd_parameter_default_values(factory, community_context, FT, cache).holling_half_saturation
        ),
        palatability_matrix = NoDefault(),
        assimilation_matrix = NoDefault(),
        optimum_predator_prey_ratio = FnDefault((factory, community_context, FT, cache) ->
            nipizd_parameter_default_values(factory, community_context, FT, cache).optimum_predator_prey_ratio
        ),
        specificity = FnDefault((factory, community_context, FT, cache) ->
            nipizd_parameter_default_values(factory, community_context, FT, cache).specificity
        ),
        protection = FnDefault((factory, community_context, FT, cache) ->
            nipizd_parameter_default_values(factory, community_context, FT, cache).protection
        ),
        assimilation_efficiency = FnDefault((factory, community_context, FT, cache) ->
            nipizd_parameter_default_values(factory, community_context, FT, cache).assimilation_efficiency
        ),
    )
end

@inline function nipizd_parameter_default_values(
    factory::NiPiZDFactory, community_context::CommunityContext, ::Type{FT}, cache::Dict{Symbol,Any}
) where {FT}
    return get!(cache, :nipizd_parameter_default_values) do
        nipizd_parameter_defaults(factory, community_context, FT)
    end
end


function nipizd_parameter_defaults(
    ::NiPiZDFactory, community_context::CommunityContext, ::Type{FT}
) where {FT}
    detritus_remineralization = FT(0.1213 / 86400)
    mortality_export_fraction = FT(0.2)

    # Group vectors (length = community_context.n_total)
    linear_mortality = fill(FT(8e-7), community_context.n_total)
    quadratic_mortality = resolve_diameter_indexed_vector(
        FT, community_context.diameters, community_context.consumer_param_indices, FT(1e-6); default=FT(0)
    )

    maximum_growth_rate = resolve_diameter_indexed_vector(
        FT,
        community_context.diameters,
        community_context.producer_param_indices,
        AllometricParam(PowerLaw(); prefactor=FT(2 / 86400), exponent=FT(-0.15));
        default=FT(0),
    )

    nutrient_half_saturation = resolve_diameter_indexed_vector(
        FT,
        community_context.diameters,
        community_context.producer_param_indices,
        AllometricParam(PowerLaw(); prefactor=FT(0.17), exponent=FT(0.27));
        default=FT(0),
    )

    alpha = resolve_diameter_indexed_vector(
        FT, community_context.diameters, community_context.producer_param_indices, FT(0.1953 / 86400); default=FT(0)
    )

    maximum_predation_rate = resolve_diameter_indexed_vector(
        FT,
        community_context.diameters,
        community_context.consumer_param_indices,
        AllometricParam(PowerLaw(); prefactor=FT(30.84 / 86400), exponent=FT(-0.16));
        default=FT(0),
    )

    holling_half_saturation = resolve_diameter_indexed_vector(
        FT, community_context.diameters, community_context.consumer_param_indices, FT(5.0); default=FT(0)
    )

    # --- Default interaction matrices (internal trait defaults) ------------

    optimum_predator_prey_ratio = resolve_diameter_indexed_vector(
        FT, community_context.diameters, community_context.consumer_param_indices, FT(10.0); default=FT(0)
    )
    specificity = resolve_diameter_indexed_vector(
        FT, community_context.diameters, community_context.consumer_param_indices, FT(0.3); default=FT(0)
    )
    protection = resolve_diameter_indexed_vector(
        FT, community_context.diameters, community_context.consumer_param_indices, FT(1.0); default=FT(0)
    )
    assimilation_efficiency = resolve_diameter_indexed_vector(
        FT, community_context.diameters, community_context.consumer_param_indices, FT(0.32); default=FT(0)
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
    ::NiPiZDFactory, community_context::CommunityContext, params::NamedTuple
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

@inline function _derive_assimilation_matrix(
    ::NiPiZDFactory, community_context::CommunityContext, params::NamedTuple
)
    FT = community_context.FT
    return assimilation_efficiency_matrix_binary_axes(
        FT;
        assimilation_efficiency=params.assimilation_efficiency,
        consumer_indices=community_context.consumer_indices,
        prey_indices=community_context.prey_indices,
    )
end

function derived_matrix_specs(::NiPiZDFactory)
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
