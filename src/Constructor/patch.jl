"""Utilities for ergonomically overriding factory configuration.

Agate uses immutable configuration objects (primarily `NamedTuple`s and lightweight
containers around `NamedTuple`s). This file provides a small, non-mutating
"patch" layer that creates overridden copies and (for `NamedTuple`s) checks that
keys exist to prevent silent typos.
"""

# -----------------------------------------------------------------------------
# NamedTuple patching
# -----------------------------------------------------------------------------

"""
    patch(nt::NamedTuple; kwargs...) -> NamedTuple

Return a copy of `nt` with overridden fields.

For fields whose current value is a `NamedTuple`, providing a `NamedTuple` override
performs a *nested* merge (and validates that nested keys exist).

Top-level keys are validated to exist in `nt` to catch typos early.
"""
function patch(nt::NamedTuple; kwargs...)
    isempty(kwargs) && return nt

    for k in keys(kwargs)
        hasproperty(nt, k) || throw(ArgumentError("patch: key $(k) is not present in this NamedTuple"))
    end

    overrides = Pair{Symbol, Any}[]
    for (k, v) in kwargs
        old = getproperty(nt, k)

        if old isa NamedTuple && v isa NamedTuple
            # Nested validation: prevent silent creation of unknown nested keys.
            for nk in keys(v)
                hasproperty(old, nk) || throw(ArgumentError("patch: key $(k).$(nk) is not present in this NamedTuple"))
            end
            push!(overrides, k => merge(old, v))
        else
            push!(overrides, k => v)
        end
    end

    return merge(nt, (; overrides...))
end

# -----------------------------------------------------------------------------
# Parameter container patching
# -----------------------------------------------------------------------------

"""
    patch(pft::PFTSpecification; kwargs...) -> PFTSpecification

Return a copy of `pft` with overridden fields.

`PFTSpecification` is intentionally flexible: new keys are allowed.
"""
patch(pft::PFTSpecification; kwargs...) = PFTSpecification(; pft.data..., kwargs...)

"""
    patch(ms::ModelSpecification; kwargs...) -> ModelSpecification

Return a copy of `ms` with overridden fields.

This is intended for lightweight experiments and tests. New keys are allowed.
"""
patch(ms::ModelSpecification; kwargs...) = ModelSpecification((; ms.data..., kwargs...))

"""    update_community(community::NamedTuple; kwargs...) -> NamedTuple

Return a copy of `community` with overridden *structural* fields.

This helper is intended for structural edits only (e.g. changing group counts or
size-class diameter ranges) without mutating the original container.

Parameter values should be overridden via `update_registry`. To remove API overlap,
`update_community` disallows overriding any `pft` field (which can hold per-PFT
parameter overrides).
"""
function update_community(community::NamedTuple; kwargs...)
    isempty(kwargs) && return community

    for k in keys(kwargs)
        hasproperty(community, k) || throw(ArgumentError("update_community: group $(k) is not present in community"))
    end

    overrides = Pair{Symbol, Any}[]
    for (k, v) in kwargs
        old = getproperty(community, k)

        # Disallow parameter overrides via PFT containers.
        if v isa PFTSpecification
            throw(ArgumentError("update_community: cannot override PFTSpecification via update_community; use update_registry for parameter overrides"))
        end

        if old isa NamedTuple && v isa NamedTuple
            # Nested validation: prevent silent creation of unknown nested keys.
            for nk in keys(v)
                hasproperty(old, nk) || throw(ArgumentError("update_community: key $(k).$(nk) is not present in this group specification"))
                nk === :pft && throw(ArgumentError("update_community: cannot override `pft` via update_community; use update_registry for parameter overrides"))
            end
            push!(overrides, k => merge(old, v))
        else
            # Reject any attempt to smuggle `pft` through a replacement NamedTuple.
            if v isa NamedTuple && hasproperty(v, :pft)
                throw(ArgumentError("update_community: cannot override `pft` via update_community; use update_registry for parameter overrides"))
            end
            push!(overrides, k => v)
        end
    end

    return merge(community, (; overrides...))
end

"""    extend_community(community::NamedTuple; kwargs...) -> NamedTuple

Return a copy of `community` with **additional** groups appended.

Unlike `update_community`, this helper allows *new* top-level keys. It errors if a
provided key already exists to prevent accidental overwrites.
"""
function extend_community(community::NamedTuple; kwargs...)
    isempty(kwargs) && return community
    for k in keys(kwargs)
        hasproperty(community, k) && throw(ArgumentError("extend_community: group $(k) already exists; use update_community to modify it"))
    end
    return merge(community, (; kwargs...))
end

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Convenience helpers
# -----------------------------------------------------------------------------

"""
    update_community(community::NamedTuple, group::Symbol; kwargs...) -> NamedTuple

Return a copy of `community` with the group entry `group` updated.

This helper is intended for *structural* edits (e.g. changing group counts or
size-class diameters) without mutating the original container.

Behavior:
- Keyword arguments must match existing keys in `community.<group>`.
- Overriding `pft` (or any per-PFT parameters) is disallowed; use `update_registry`.
"""
function update_community(community::NamedTuple, group::Symbol; kwargs...)
    hasproperty(community, group) || throw(ArgumentError("update_community: group $(group) is not present in community"))
    spec = getproperty(community, group)
    spec isa NamedTuple || throw(ArgumentError("update_community: community.$(group) must be a NamedTuple, got $(typeof(spec))"))

    isempty(kwargs) && return community

    for k in keys(kwargs)
        k === :pft && throw(ArgumentError("update_community: cannot override `pft` via update_community; use update_registry for parameter overrides"))
        hasproperty(spec, k) || throw(ArgumentError("update_community: key $(k) is not present in community.$(group); use update_registry for parameter overrides"))
    end

    patched_spec = patch(spec; kwargs...)

    return merge(community, (; (group => patched_spec),))
end

"""
    update_dynamics(dynamics::NamedTuple; kwargs...) -> NamedTuple

Return a copy of `dynamics` with overridden entries.

Keys are validated to exist to prevent silent typos.
"""
function update_dynamics(dynamics::NamedTuple; kwargs...)
    isempty(kwargs) && return dynamics
    for k in keys(kwargs)
        hasproperty(dynamics, k) || throw(ArgumentError("update_dynamics: key $(k) is not present in dynamics"))
    end
    return merge(dynamics, (; kwargs...))
end

"""    extend_dynamics(dynamics::NamedTuple; kwargs...) -> NamedTuple

Return a copy of `dynamics` with **new** entries appended.

Errors if a provided key already exists to prevent accidental overwrites.
"""
function extend_dynamics(dynamics::NamedTuple; kwargs...)
    isempty(kwargs) && return dynamics
    for k in keys(kwargs)
        hasproperty(dynamics, k) && throw(ArgumentError("extend_dynamics: key $(k) already exists; use update_dynamics to modify it"))
    end
    return merge(dynamics, (; kwargs...))
end
