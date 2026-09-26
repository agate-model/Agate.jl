export default_components
export default_processes
export definition_version
export plankton_roles
export ModelSetting
export setting_definitions

"""One scientific model-family setting that is resolved at construction but not bound to a process slot."""
struct ModelSetting{Default}
    default::Default
    domain::Symbol

    function ModelSetting(default::Default, domain::Symbol) where {Default}
        domain in (:any, :finite, :nonnegative, :positive, :unit_interval) || throw(
            ArgumentError(
                "ModelSetting domain must be :any, :finite, :nonnegative, :positive, or :unit_interval"
            ),
        )
        return new{Default}(default, domain)
    end
end

ModelSetting(default; domain::Symbol=:finite) = ModelSetting(default, domain)

"""Scientific settings for a named model family that are not process-bound parameters."""
setting_definitions(::AbstractModelFamily) = (;)

"""Map user-facing plankton roles to logical components, e.g. `phytoplankton => :P`."""
plankton_roles(::AbstractModelFamily) = throw(
    ArgumentError("No method `plankton_roles(family)` is defined for this model family.")
)

"""Canonical logical components for a named model family.

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

"""Canonical named scientific processes for a named model family.

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

"""Scientific definition version for a named model family.

Bump this version whenever family science, defaults, derivation algorithms, or canonical
definition structure changes in a way that should invalidate durable recipe replay.
"""
function definition_version(::AbstractModelFamily)::VersionNumber
    throw(ArgumentError("No method `definition_version(family)` is defined for this model family."))
end
