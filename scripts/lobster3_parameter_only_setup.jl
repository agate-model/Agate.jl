# Parameter-only FrankenLOBSTER approximation of the supplied LOBSTER3 model.
#
# FrankenLOBSTER defaults already match the LOBSTER3 3-D allometries for:
#   P maximum growth, P nitrate affinity, H maximum uptake, H DOM affinity,
#   P/Z/H mortality, Z/H assimilation, C:N, and chlorophyll:N.
# Only deliberate departures from those defaults are specified below.
#
# Temperature experiments:
#   :temperature_off => Q10(P1, P2) = (1.0, 1.0)
#   :p1_temperature  => Q10(P1, P2) = (1.88, 1.0)
#
# Remaining functional-form approximations:
# - LOBSTER3 mass-action grazing is approached with K_G >> prey and
#   g_max(d) = K_G * 15.9 * V(d)^(-0.16) / day.
# - The inactive LOBSTER3 NH4 growth branch is made negligible with a very
#   large NH4 half-saturation.
# - LOBSTER3 Monod light (K=55 W m^-2) is matched at its 50% point by the
#   FrankenLOBSTER exponential-saturation light response.

using Agate
using OceanBioME
using Oceananigans
using Oceananigans.Units: day
using Agate.Library.Allometry: AllometricParam, ConstantParam, PowerLaw
using OceanBioME.Models.NutrientsPlanktonDetritusModels: DissolvedParticulate, LOBSTER

const FrankenLOBSTER = Agate.Models.FrankenLOBSTER

const GRAZING_K = 100.0       # >99% of mass-action rate for total prey < 1 mmol N m^-3
const NH4_OFF_K = 1.0e12
const LIGHT_K = 55 / log(2)   # exponential response = 0.5 at PAR = 55 W m^-2

function temperature_q10(experiment::Symbol)
    experiment === :temperature_off && return (P_1=1.0, P_2=1.0)
    experiment === :p1_temperature && return (P_1=1.88, P_2=1.0)
    throw(ArgumentError("temperature_experiment must be :temperature_off or :p1_temperature"))
end

const BASE_PARAMETERS = (
    # LOBSTER3 currently uses nitrate-supported P growth only.
    ammonia_half_saturation=ConstantParam(NH4_OFF_K),
    nitrate_ammonia_inhibition=0.0,

    # LOBSTER3 uses Monod(PAR; K=55); FrankenLOBSTER uses exponential saturation.
    light_half_saturation=(P_1=LIGHT_K, P_2=LIGHT_K),

    # Disabled in the supplied LOBSTER3 setup.
    phytoplankton_exudation_fraction=(P_1=0.0, P_2=0.0),
    zooplankton_excretion_rate=(Z_1=0.0, Z_2=0.0),

    # Z1 -> P1 + H1; Z2 -> P2. Scale g_max with K_G so the Holling response
    # approaches LOBSTER3's mass-action g(d) * prey * predator formulation.
    maximum_predation_rate=AllometricParam(
        PowerLaw(); prefactor=GRAZING_K * 15.9 / day, exponent=-0.16
    ),
    grazing_half_saturation=(Z_1=GRAZING_K, Z_2=GRAZING_K),
    palatability_matrix=[1.0 0.0 1.0; 0.0 1.0 0.0],
)

lobster3_parameters(; temperature_experiment=:temperature_off) = merge(
    BASE_PARAMETERS,
    (temperature_q10=temperature_q10(temperature_experiment),),
)

function lobster3_like_bgc(grid; temperature_experiment=:temperature_off)
    plankton = FrankenLOBSTER.construct(;
        grid,
        parameters=lobster3_parameters(; temperature_experiment),
    )

    # The only LOBSTER detritus default changed by the supplied LOBSTER3 setup.
    detritus = DissolvedParticulate(
        grid;
        dissolved_remineralisation_rate=0.0,
    )

    return LOBSTER(grid; plankton, detritus)
end

# Directly runnable one-cell setups for the two intended experiments.
grid = RectilinearGrid(CPU(); size=(1, 1, 1), extent=(1, 1, 1))
bgc_temperature_off = lobster3_like_bgc(grid; temperature_experiment=:temperature_off)
bgc_p1_temperature = lobster3_like_bgc(grid; temperature_experiment=:p1_temperature)
