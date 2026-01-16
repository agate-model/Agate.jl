"""Parameter registry for the DARWIN model.

Defaults in this file are the **single source of truth** for DARWIN parameters.
Example/spec containers may provide structural information and user overrides only.
"""

using Agate.Parameters: ParamRegistry
using Agate.Parameters: scalar_param, vector_param, matrix_param
using Agate.Parameters: GroupVec
using Agate.Models: default_palatability_provider, default_assimilation_provider
using Agate.Library.Allometry: AllometricParam, PowerLaw

import Agate.Parameters: parameter_registry

function parameter_registry(::DarwinFactory)
    groups = (:Z, :P)
    return ParamRegistry([
        # --- Stoichiometry + partitioning scalars --------------------------
        scalar_param(:nitrogen_to_carbon, "Nitrogen:carbon ratio.", 0.15; missing_policy=:fail),
        scalar_param(:phosphorus_to_carbon, "Phosphorus:carbon ratio.", 0.009; missing_policy=:fail),
        scalar_param(:DOM_POM_fractionation, "Fractionation between DOM and POM.", 0.45; missing_policy=:fail),

        # Remineralization rates
        scalar_param(:POC_remineralization, "POC remineralization rate.", 0.1213 / 86400; missing_policy=:fail),
        scalar_param(:DOC_remineralization, "DOC remineralization rate.", 0.1213 / 86400; missing_policy=:fail),
        scalar_param(:PON_remineralization, "PON remineralization rate.", 0.1213 / 86400; missing_policy=:fail),
        scalar_param(:DON_remineralization, "DON remineralization rate.", 0.1213 / 86400; missing_policy=:fail),
        scalar_param(:POP_remineralization, "POP remineralization rate.", 0.1213 / 86400; missing_policy=:fail),
        scalar_param(:DOP_remineralization, "DOP remineralization rate.", 0.1213 / 86400; missing_policy=:fail),

        # --- Plankton vectors ----------------------------------------------
        vector_param(
            :linear_mortality,
            "Linear mortality rate.",
            GroupVec(groups; Z=8e-7, P=8e-7);
            missing_policy=:fail,
        ),
        vector_param(
            :quadratic_mortality,
            "Quadratic mortality rate.",
            GroupVec(groups; Z=1e-6, P=0.0);
            missing_policy=:fail,
        ),

        # Photosynthesis / nutrient limitation
        vector_param(
            :maximum_growth_rate,
            "Maximum phytoplankton growth rate.",
            GroupVec(groups; P=AllometricParam(PowerLaw(); prefactor=2 / 86400, exponent=-0.15), Z=0.0);
            missing_policy=:fail,
        ),
        vector_param(
            :half_saturation_DIN,
            "DIN half-saturation.",
            GroupVec(groups; P=AllometricParam(PowerLaw(); prefactor=0.17, exponent=0.27), Z=0.0);
            missing_policy=:fail,
        ),
        vector_param(
            :half_saturation_PO4,
            "PO4 half-saturation.",
            GroupVec(groups; P=AllometricParam(PowerLaw(); prefactor=0.17, exponent=0.27), Z=0.0);
            missing_policy=:fail,
        ),
        vector_param(:alpha, "Smith-style light-limitation parameter.", 0.0 / 86400; missing_policy=:fail),
        vector_param(
            :photosynthetic_slope,
            "Photosynthetic slope.",
            GroupVec(groups; P=0.46e-5, Z=0.0);
            missing_policy=:fail,
        ),
        vector_param(
            :chlorophyll_to_carbon_ratio,
            "Chlorophyll:carbon ratio.",
            GroupVec(groups; P=0.1, Z=0.0);
            missing_policy=:fail,
        ),

        # --- Grazing community vectors -------------------------------------
        vector_param(
            :maximum_predation_rate,
            "Maximum grazing rate.",
            GroupVec(groups; Z=AllometricParam(PowerLaw(); prefactor=30.84 / 86400, exponent=-0.16), P=0.0);
            missing_policy=:fail,
        ),
        vector_param(
            :holling_half_saturation,
            "Holling half-saturation constant.",
            GroupVec(groups; Z=5.0, P=0.0);
            missing_policy=:fail,
        ),

        # --- Traits for matrix defaults ------------------------------------
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

        # --- Interaction matrices ------------------------------------------
        matrix_param(:palatability_matrix, "Predator-prey palatability matrix.", default_palatability_provider(); missing_policy=:fail),
        matrix_param(:assimilation_matrix, "Predator-prey assimilation-efficiency matrix.", default_assimilation_provider(); missing_policy=:fail),
    ])
end
