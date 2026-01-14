"""Parameter registry for the NiPiZD model.

Defaults in this file are the **single source of truth** for NiPiZD parameters.
Example/spec containers (e.g. `default_parameter_args`) should not duplicate defaults;
they may only provide structural info (sizes/diameters) and user overrides.
"""

using ...Parameters: ParamSpec, ParamRegistry
using ...Parameters: default_palatability_matrix, default_assimilation_efficiency_matrix

using ...Library.Allometry: AllometricParam, PowerLaw

import ...Parameters: parameter_registry

function parameter_registry(::NiPiZDFactory)
    return ParamRegistry([
        # --- BGC scalars -----------------------------------------------------
        ParamSpec(:detritus_remineralization, "Detritus remineralization rate.",
                 0.1213 / 86400; scope=:fail),
        ParamSpec(:mortality_export_fraction, "Fraction of mortality exported (rest becomes detritus).",
                 0.2; scope=:fail),

        # --- Plankton vectors (rates, half-saturations) ----------------------
        ParamSpec(:linear_mortality, "Linear mortality rate (applies to all plankton unless overridden).",
                 8e-7 ; scope=:fail),
        ParamSpec(:quadratic_mortality, "Quadratic mortality rate (defaults only for zooplankton).",
                 (Z = 1e-6 ,); scope=:zero_silent),
        ParamSpec(:maximum_growth_rate, "Maximum phytoplankton growth rate (allometric for phyto).",
                 (P = AllometricParam(PowerLaw(); prefactor = 2 / 86400, exponent = -0.15),); scope=:zero_silent),
        ParamSpec(:nutrient_half_saturation, "Half-saturation constant for nutrient limitation (allometric for phyto).",
                 (P = AllometricParam(PowerLaw(); prefactor = 0.17, exponent = 0.27),); scope=:zero_silent),
        ParamSpec(:alpha, "Smith-style light-limitation parameter (phyto only).",
                 (P = 0.1953 / 86400,); scope=:zero_silent),

        # --- Grazing community vectors --------------------------------------
        ParamSpec(:maximum_predation_rate, "Maximum grazing rate (allometric for zoo; phyto entries default to 0).",
                 (Z = AllometricParam(PowerLaw(); prefactor = 30.84 / 86400, exponent = -0.16),); scope=:zero_silent),
        ParamSpec(:holling_half_saturation, "Holling type-II half-saturation constant (zoo only).",
                 (Z = 5.0,); scope=:zero_silent),

        # --- Traits used to build interaction matrices (computed on CPU) ----
        ParamSpec(:can_eat, "Binary predator flag used to build palatability / assimilation matrices.",
                 (Z = true, P = false); scope=:fail, kind=:bool),
        ParamSpec(:can_be_eaten, "Binary prey flag used to build assimilation matrix.",
                 (Z = false, P = true); scope=:fail, kind=:bool),
        ParamSpec(:optimum_predator_prey_ratio, "Preferred predator:prey size ratio (zoo only).",
                 (Z = 10.0,); scope=:zero_silent),
        ParamSpec(:specificity, "Palatability specificity (zoo only).",
                 (Z = 0.3,); scope=:zero_silent),
        ParamSpec(:protection, "Prey protection factor (zoo only; prey entries default to 0).",
                 (Z = 1.0,); scope=:zero_silent),
        ParamSpec(:assimilation_efficiency, "Assimilation efficiency used to build assimilation matrix (zoo only).",
                 (Z = 0.32,); scope=:zero_silent),

        # --- Interaction matrices -------------------------------------------
        ParamSpec(:palatability_matrix, "Predator–prey palatability matrix.",
                 default_palatability_matrix; scope=:fail),
        ParamSpec(:assimilation_efficiency_matrix, "Predator–prey assimilation efficiency matrix.",
                 default_assimilation_efficiency_matrix; scope=:fail),
    ])
end

# -----------------------------------------------------------------------------
# Parameter placeholders for equation authoring
# -----------------------------------------------------------------------------

# Define and export `ParamVar` placeholders with the same names as all registry
# entries. These are used only at construction time when building symbolic
# `Equation`s; runtime kernels see only numeric arrays/scalars.
using ...Library.Equations: declare_parameter_vars!

const _nipizd_param_names = map(s -> s.name, parameter_registry(NiPiZDFactory()).specs)

# Declare parameter placeholders both in the parent model module and in the
# `Tracers` submodule where equations are authored.
declare_parameter_vars!(@__MODULE__, _nipizd_param_names; export_vars=true)
if isdefined(@__MODULE__, :Tracers)
    declare_parameter_vars!(Tracers, _nipizd_param_names; export_vars=false)
end
