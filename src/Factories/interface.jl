export default_components
export default_processes

"""Canonical logical components for a registered model family.

Returns a named collection whose keys are stable model component identities and
whose values describe intrinsic component structure.
"""
function default_components(::AbstractBGCFactory)
    throw(
        ArgumentError(
            "No method `default_components(factory)` is defined for this model family."
        ),
    )
end

"""Canonical named scientific processes for a registered model family.

The keys are stable process-instance identities. Process declarations describe
scientific topology and are normalized before runtime realization.
"""
function default_processes(::AbstractBGCFactory)
    throw(
        ArgumentError(
            "No method `default_processes(factory)` is defined for this model family."
        ),
    )
end
