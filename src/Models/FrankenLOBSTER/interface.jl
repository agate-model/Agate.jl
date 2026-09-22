using Adapt: adapt
import Adapt: adapt_structure

import Oceananigans.Biogeochemistry:
    biogeochemical_drift_velocity,
    required_biogeochemical_auxiliary_fields,
    required_biogeochemical_tracers

"""Internal OceanBioME plankton component backed by a compiled Agate realization.

`OwnedTracers` are the living P/Z/B tracers registered by the surrounding LOBSTER model.
`ExternalTracers` are OceanBioME-owned resource fields read by compiled Agate processes.
Both are encoded in the type so architecture adaptation does not carry Symbol metadata into
runtime storage.
"""
struct FrankenLOBSTERPlankton{Runtime,OwnedTracers,ExternalTracers}
    runtime::Runtime
end

FrankenLOBSTERPlankton(runtime, owned::Tuple, external::Tuple) =
    FrankenLOBSTERPlankton{typeof(runtime),owned,external}(runtime)

@inline required_biogeochemical_tracers(
    ::FrankenLOBSTERPlankton{Runtime,OwnedTracers}
) where {Runtime,OwnedTracers} = OwnedTracers

@inline required_biogeochemical_auxiliary_fields(plankton::FrankenLOBSTERPlankton) =
    required_biogeochemical_auxiliary_fields(plankton.runtime)

@inline biogeochemical_drift_velocity(plankton::FrankenLOBSTERPlankton, tracer::Val) =
    biogeochemical_drift_velocity(plankton.runtime, tracer)

@inline (plankton::FrankenLOBSTERPlankton)(tracer::Val, args...) =
    plankton.runtime(tracer, args...)

@inline external_tracers(
    ::FrankenLOBSTERPlankton{Runtime,OwnedTracers,ExternalTracers}
) where {Runtime,OwnedTracers,ExternalTracers} = ExternalTracers

@inline function adapt_structure(
    to,
    plankton::FrankenLOBSTERPlankton{Runtime,OwnedTracers,ExternalTracers},
) where {Runtime,OwnedTracers,ExternalTracers}
    return FrankenLOBSTERPlankton(
        adapt(to, plankton.runtime), OwnedTracers, ExternalTracers
    )
end
