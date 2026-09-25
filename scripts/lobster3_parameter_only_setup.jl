using Agate
using OceanBioME
using Oceananigans
using Oceananigans.Units: day
using Agate.Library.Allometry: AllometricParam, ConstantParam, PowerLaw
using OceanBioME.Models.NutrientsPlanktonDetritusModels: DissolvedParticulate, LOBSTER

const FrankenLOBSTER = Agate.Models.FrankenLOBSTER

# to make grazing linear in prey for total prey < 1 mmol N m^-3, we need to scale g_max with K_G
const GRAZING_K = 100.0       # >99% of mass-action rate for total prey < 1 mmol N m^-3
const NH4_OFF_K = 1.0e12

# to convert from LOBSTER3's Monod(PAR; K=55) to FrankenLOBSTER's exponential saturation
const LIGHT_K = 55 / log(2)   # exponential response = 0.5 at PAR = 55 W m^-2

# temperature response on/off flag (1 = exponential of 0 for Q10, e.g. no temperature response )
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
    # approaches LOBSTER3's g(d) * prey * predator formulation.
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
    # Agate model passed to OceanBioME LOBSTER
    plankton = FrankenLOBSTER.construct(;
        grid,
        parameters=lobster3_parameters(; temperature_experiment),
    )

    # There is no DOM remin in LOBSTER3 so this has to be redefined on the LOBSTER side:
    detritus = DissolvedParticulate(
        grid;
        dissolved_remineralisation_rate=0.0,
    )

    return LOBSTER(grid; plankton, detritus)
end

# Dummy example of model run:
grid = RectilinearGrid(CPU(); size=(1, 1, 1), extent=(1, 1, 1))
# with temperature off for both P1 and P2
bgc_temperature_off = lobster3_like_bgc(grid; temperature_experiment=:temperature_off)
# with temperature on for P1 and off for P2
bgc_p1_temperature = lobster3_like_bgc(grid; temperature_experiment=:p1_temperature)
