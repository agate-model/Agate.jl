using ...ModelFamilies: AbstractModelFamily
using ...Components: Plankton, Pool

import ...ModelFamilies: default_components, default_processes, definition_version
import ...Construction: family_id, registered_family

"""Registered family for the Agate-generated FrankenLOBSTER living community."""
struct FrankenLOBSTERFamily <: AbstractModelFamily end

family_id(::FrankenLOBSTERFamily) = :FrankenLOBSTER
registered_family(::Val{:FrankenLOBSTER}) = FrankenLOBSTERFamily()
definition_version(::FrankenLOBSTERFamily)::VersionNumber = v"0.1.0"

"""LOBSTER3-like default living-community size structure."""
const DEFAULT_SIZE_STRUCTURE = (
    phytoplankton=(P=(n=2, min_esd=0.6, max_esd=1.2, spacing=:linear),),
    zooplankton=(Z=(n=2, min_esd=6.0, max_esd=12.0, spacing=:linear),),
    bacterioplankton=(B=(n=1, min_esd=0.6, max_esd=0.6, spacing=:linear),),
)

# NO3, NH4, and DOM are OceanBioME-owned state in FrankenLOBSTER. They are represented
# here so Agate processes can use the same named resource identities when compiling the
# living-community equations. The plankton adapter exposes only P/Z/B as owned tracers.
const FRANKENLOBSTER_COMPONENTS = (
    NO₃=Pool(:nitrogen),
    NH₄=Pool(:nitrogen),
    DOM=Pool(:nitrogen),
    P=Plankton(;
        states=(nitrogen=:nitrogen,),
        reference_state=:nitrogen,
        size_structure=DEFAULT_SIZE_STRUCTURE.phytoplankton.P,
    ),
    Z=Plankton(;
        states=(nitrogen=:nitrogen,),
        reference_state=:nitrogen,
        size_structure=DEFAULT_SIZE_STRUCTURE.zooplankton.Z,
    ),
    B=Plankton(;
        states=(nitrogen=:nitrogen,),
        reference_state=:nitrogen,
        size_structure=DEFAULT_SIZE_STRUCTURE.bacterioplankton.B,
    ),
)

"""Canonical logical components for FrankenLOBSTER."""
default_components(::FrankenLOBSTERFamily) = FRANKENLOBSTER_COMPONENTS

# Cycle 1 establishes the construction/integration boundary. Scientific process definitions
# are added in the following cycles without changing the adapter ownership contract.
const FRANKENLOBSTER_PROCESSES = (;)

"""Canonical named scientific processes for FrankenLOBSTER."""
default_processes(::FrankenLOBSTERFamily) = FRANKENLOBSTER_PROCESSES
