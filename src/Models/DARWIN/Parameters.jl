"""Parameter registry for the DARWIN model.

Defaults in this file are the **single source of truth** for DARWIN parameters.
Example/spec containers may provide structural information and user overrides only.
"""

using ...Parameters: ParamSpec, ParamRegistry
using ...Parameters: default_palatability_matrix, default_assimilation_matrix

using ...Library.Allometry: AllometricParam, PowerLaw

import ...Parameters: parameter_registry

function parameter_registry(::DarwinFactory)
    return ParamRegistry([
        # --- Stoichiometry + partitioning scalars --------------------------
        ParamSpec(:nitrogen_to_carbon, "Nitrogen:carbon ratio.",
                 0.15; missing_policy=:fail),
        ParamSpec(:phosphorus_to_carbon, "Phosphorus:carbon ratio.",
                 0.009; missing_policy=:fail),
        ParamSpec(:DOM_POM_fractionation, "Fractionation between DOM and POM.",
                 0.45; missing_policy=:fail),

        # Remineralization rates
        ParamSpec(:POC_remineralization, "POC remineralization rate.",
                 0.1213 / 86400; missing_policy=:fail),
        ParamSpec(:DOC_remineralization, "DOC remineralization rate.",
                 0.1213 / 86400; missing_policy=:fail),
        ParamSpec(:PON_remineralization, "PON remineralization rate.",
                 0.1213 / 86400; missing_policy=:fail),
        ParamSpec(:DON_remineralization, "DON remineralization rate.",
                 0.1213 / 86400; missing_policy=:fail),
        ParamSpec(:POP_remineralization, "POP remineralization rate.",
                 0.1213 / 86400; missing_policy=:fail),
        ParamSpec(:DOP_remineralization, "DOP remineralization rate.",
                 0.1213 / 86400; missing_policy=:fail),

        # --- Plankton vectors ----------------------------------------------
        ParamSpec(:linear_mortality, "Linear mortality rate.",
                 8e-7 ; missing_policy=:fail),
        ParamSpec(:quadratic_mortality, "Quadratic mortality rate (zoo only).",
                 (Z = 1e-6 ,); missing_policy=:zero_silent),

        # Photosynthesis / nutrient limitation
        ParamSpec(:maximum_growth_rate, "Maximum phytoplankton growth rate (allometric for phyto).",
                 (P = AllometricParam(PowerLaw(); prefactor = 2 / 86400, exponent = -0.15),); missing_policy=:zero_silent),
        ParamSpec(:half_saturation_DIN, "DIN half-saturation (allometric for phyto).",
                 (P = AllometricParam(PowerLaw(); prefactor = 0.17, exponent = 0.27),); missing_policy=:zero_silent),
        ParamSpec(:half_saturation_PO4, "PO4 half-saturation (allometric for phyto).",
                 (P = AllometricParam(PowerLaw(); prefactor = 0.17, exponent = 0.27),); missing_policy=:zero_silent),
        ParamSpec(:alpha, "Smith-style light-limitation parameter.",
                 0.0 / 86400; missing_policy=:fail),
        ParamSpec(:photosynthetic_slope, "Photosynthetic slope (phyto only).",
                 (P = 0.46e-5,); missing_policy=:zero_silent),
        ParamSpec(:chlorophyll_to_carbon_ratio, "Chlorophyll:carbon ratio (phyto only).",
                 (P = 0.1,); missing_policy=:zero_silent),

        # --- Grazing community vectors -------------------------------------
        ParamSpec(:maximum_predation_rate, "Maximum grazing rate (allometric for zoo).",
                 (Z = AllometricParam(PowerLaw(); prefactor = 30.84 / 86400, exponent = -0.16),); missing_policy=:zero_silent),
        ParamSpec(:holling_half_saturation, "Holling half-saturation constant (zoo only).",
                 (Z = 5.0,); missing_policy=:zero_silent),

        # --- Traits for matrix defaults ------------------------------------
        ParamSpec(:can_eat, "Binary predator flag used for default interactions.",
                 (Z = true, P = false); missing_policy=:fail, value_kind=:bool),
        ParamSpec(:can_be_eaten, "Binary prey flag used for default interactions.",
                 (Z = false, P = true); missing_policy=:fail, value_kind=:bool),
        ParamSpec(:optimum_predator_prey_ratio, "Preferred predator:prey size ratio.",
                 (Z = 10.0,); missing_policy=:zero_silent),
        ParamSpec(:specificity, "Palatability specificity.",
                 (Z = 0.3,); missing_policy=:zero_silent),
        ParamSpec(:protection, "Prey protection factor.",
                 (Z = 1.0,); missing_policy=:zero_silent),
        ParamSpec(:assimilation_efficiency, "Assimilation efficiency used to build the default assimilation-efficiency matrix.",
                 (Z = 0.32,); missing_policy=:zero_silent),

        # --- Interaction matrices ------------------------------------------
        ParamSpec(:palatability_matrix, "Predator–prey palatability matrix.",
                 default_palatability_matrix; missing_policy=:fail),
        ParamSpec(:assimilation_matrix, "Predator–prey assimilation matrix.",
                 default_assimilation_matrix; missing_policy=:fail),
    ])
end
