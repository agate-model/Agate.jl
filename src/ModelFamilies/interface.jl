export default_components
export default_processes
export definition_version

"""Canonical logical components for a registered model family.

Returns a named collection whose keys are stable model component identities and
whose values describe intrinsic component structure.
"""
function default_components(::AbstractModelFamily)
    throw(
        ArgumentError(
            "No method `default_components(family)` is defined for this model family."
        ),
    )
end

"""Canonical named scientific processes for a registered model family.

The keys are stable process-instance identities. Process declarations describe
scientific topology and are canonicalized before runtime realization.
"""
function default_processes(::AbstractModelFamily)
    throw(
        ArgumentError(
            "No method `default_processes(family)` is defined for this model family."
        ),
    )
end

"""Scientific definition version for a registered model family.

Bump this version whenever family science, defaults, derivation algorithms, or canonical
definition structure changes in a way that should invalidate durable recipe replay.
"""
function definition_version(::AbstractModelFamily)::VersionNumber
    throw(ArgumentError("No method `definition_version(family)` is defined for this model family."))
end
