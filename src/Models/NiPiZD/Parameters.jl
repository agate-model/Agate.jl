"""Parameter registry for the NiPiZD model.

Defaults in this file are the **single source of truth** for NiPiZD parameters.
Example/spec containers (e.g. `default_community`) should not duplicate defaults;
they may only provide structural information (sizes/diameters) and user overrides.
"""

using Agate.Parameters: ParamRegistry
using Agate.Parameters: scalar_param, vector_param, matrix_param
using Agate.Parameters: GroupVec
using Agate.Models: default_palatability_provider, default_assimilation_provider
using Agate.Library.Allometry: AllometricParam, PowerLaw

import Agate.Parameters: parameter_registry

function parameter_registry(::NiPiZDFactory)
    groups = (:Z, :P)
    return ParamRegistry([
        # --- BGC scalars -----------------------------------------------------
        scalar_param(:detritus_remineralization, "Detritus remineralization rate.", 0.1213 / 86400; missing_policy=:fail),
        scalar_param(:mortality_export_fraction, "Fraction of mortality exported (rest becomes detritus).", 0.2; missing_policy=:fail),

        # --- Plankton vectors (rates, half-saturations) ----------------------
        vector_param(
            :linear_mortality,
            "Linear mortality rate (applies to all plankton unless overridden).",
            GroupVec(groups; Z=8e-7, P=8e-7);
            missing_policy=:fail,
        ),
        vector_param(
            :quadratic_mortality,
            "Quadratic mortality rate.",
            GroupVec(groups; Z=1e-6, P=0.0);
            missing_policy=:fail,
        ),
        vector_param(
            :maximum_growth_rate,
            "Maximum phytoplankton growth rate.",
            GroupVec(groups; P=AllometricParam(PowerLaw(); prefactor=2 / 86400, exponent=-0.15), Z=0.0);
            missing_policy=:fail,
        ),
        vector_param(
            :nutrient_half_saturation,
            "Half-saturation constant for nutrient limitation.",
            GroupVec(groups; P=AllometricParam(PowerLaw(); prefactor=0.17, exponent=0.27), Z=0.0);
            missing_policy=:fail,
        ),
        vector_param(
            :alpha,
            "Smith-style light-limitation parameter.",
            GroupVec(groups; P=0.1953 / 86400, Z=0.0);
            missing_policy=:fail,
        ),

        # --- Grazing community vectors --------------------------------------
        vector_param(
            :maximum_predation_rate,
            "Maximum grazing rate.",
            GroupVec(groups; Z=AllometricParam(PowerLaw(); prefactor=30.84 / 86400, exponent=-0.16), P=0.0);
            missing_policy=:fail,
        ),
        vector_param(
            :holling_half_saturation,
            "Holling type-II half-saturation constant.",
            GroupVec(groups; Z=5.0, P=0.0);
            missing_policy=:fail,
        ),

        # --- Traits used to build interaction matrices -----------------------
        vector_param(
            :can_eat,
            "Binary predator flag used to build default interaction matrices.",
            GroupVec(groups; Z=true, P=false);
            missing_policy=:fail,
            value_kind=:bool,
        ),
        vector_param(
            :can_be_eaten,
            "Binary prey flag used to build the default assimilation-efficiency matrix.",
            GroupVec(groups; Z=false, P=true);
            missing_policy=:fail,
            value_kind=:bool,
        ),
        vector_param(
            :optimum_predator_prey_ratio,
            "Preferred predator:prey size ratio.",
            GroupVec(groups; Z=10.0, P=0.0);
            missing_policy=:fail,
        ),
        vector_param(
            :specificity,
            "Palatability specificity.",
            GroupVec(groups; Z=0.3, P=0.0);
            missing_policy=:fail,
        ),
        vector_param(
            :protection,
            "Prey protection factor.",
            GroupVec(groups; Z=1.0, P=0.0);
            missing_policy=:fail,
        ),
        vector_param(
            :assimilation_efficiency,
            "Assimilation efficiency used to build the default assimilation-efficiency matrix.",
            GroupVec(groups; Z=0.32, P=0.0);
            missing_policy=:fail,
        ),

        # --- Interaction matrices -------------------------------------------
        matrix_param(:palatability_matrix, "Predator-prey palatability matrix.", default_palatability_provider(); missing_policy=:fail),
        matrix_param(:assimilation_matrix, "Predator-prey assimilation-efficiency matrix.", default_assimilation_provider(); missing_policy=:fail),
    ])
end
