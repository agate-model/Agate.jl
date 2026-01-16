"""Parameter registry for the NiPiZD model.

Defaults in this file are the **single source of truth** for NiPiZD parameters.
Example/spec containers (e.g. `default_community`) should not duplicate defaults;
they may only provide structural information (sizes/diameters) and user overrides.
"""

using Agate.Parameters: ParamRegistry
using Agate.Parameters: scalar_param, vector_param, matrix_param
using Agate.Parameters: default_palatability_provider, default_assimilation_provider
using Agate.Library.Allometry: AllometricParam, PowerLaw

import Agate.Parameters: parameter_registry

function parameter_registry(::NiPiZDFactory)
    return ParamRegistry([
        # --- BGC scalars -----------------------------------------------------
        scalar_param(:detritus_remineralization, "Detritus remineralization rate.", 0.1213 / 86400; missing_policy=:fail),
        scalar_param(:mortality_export_fraction, "Fraction of mortality exported (rest becomes detritus).", 0.2; missing_policy=:fail),

        # --- Plankton vectors (rates, half-saturations) ----------------------
        vector_param(:linear_mortality, "Linear mortality rate (applies to all plankton unless overridden).", 8e-7; missing_policy=:fail),
        vector_param(:quadratic_mortality, "Quadratic mortality rate (defaults only for zooplankton).", (Z=1e-6,); missing_policy=:zero_silent),
        vector_param(
            :maximum_growth_rate,
            "Maximum phytoplankton growth rate (allometric for phyto).",
            (P=AllometricParam(PowerLaw(); prefactor=2 / 86400, exponent=-0.15),);
            missing_policy=:zero_silent,
        ),
        vector_param(
            :nutrient_half_saturation,
            "Half-saturation constant for nutrient limitation (allometric for phyto).",
            (P=AllometricParam(PowerLaw(); prefactor=0.17, exponent=0.27),);
            missing_policy=:zero_silent,
        ),
        vector_param(:alpha, "Smith-style light-limitation parameter (phyto only).", (P=0.1953 / 86400,); missing_policy=:zero_silent),

        # --- Grazing community vectors --------------------------------------
        vector_param(
            :maximum_predation_rate,
            "Maximum grazing rate (allometric for zoo; phyto entries default to 0).",
            (Z=AllometricParam(PowerLaw(); prefactor=30.84 / 86400, exponent=-0.16),);
            missing_policy=:zero_silent,
        ),
        vector_param(:holling_half_saturation, "Holling type-II half-saturation constant (zoo only).", (Z=5.0,); missing_policy=:zero_silent),

        # --- Traits used to build interaction matrices -----------------------
        vector_param(
            :can_eat,
            "Binary predator flag used to build default interaction matrices.",
            (Z=true, P=false);
            missing_policy=:fail,
            value_kind=:bool,
        ),
        vector_param(
            :can_be_eaten,
            "Binary prey flag used to build the default assimilation-efficiency matrix.",
            (Z=false, P=true);
            missing_policy=:fail,
            value_kind=:bool,
        ),
        vector_param(:optimum_predator_prey_ratio, "Preferred predator:prey size ratio (zoo only).", (Z=10.0,); missing_policy=:zero_silent),
        vector_param(:specificity, "Palatability specificity (zoo only).", (Z=0.3,); missing_policy=:zero_silent),
        vector_param(:protection, "Prey protection factor (zoo only; prey entries default to 0).", (Z=1.0,); missing_policy=:zero_silent),
        vector_param(
            :assimilation_efficiency,
            "Assimilation efficiency used to build the default assimilation-efficiency matrix (zoo only).",
            (Z=0.32,);
            missing_policy=:zero_silent,
        ),

        # --- Interaction matrices -------------------------------------------
        matrix_param(:palatability_matrix, "Predator-prey palatability matrix.", default_palatability_provider(); missing_policy=:fail),
        matrix_param(:assimilation_matrix, "Predator-prey assimilation-efficiency matrix.", default_assimilation_provider(); missing_policy=:fail),
    ])
end
