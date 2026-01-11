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
    patch(spec::BiogeochemistrySpecification; kwargs...) -> BiogeochemistrySpecification

Return a copy of `spec` with overridden fields.

`BiogeochemistrySpecification` is intentionally flexible: new keys are allowed.
"""
patch(spec::BiogeochemistrySpecification; kwargs...) = BiogeochemistrySpecification(; spec.data..., kwargs...)

"""
    patch(ms::ModelSpecification; kwargs...) -> ModelSpecification

Return a copy of `ms` with overridden fields.

This is intended for lightweight experiments and tests. New keys are allowed.
"""
patch(ms::ModelSpecification; kwargs...) = ModelSpecification((; ms.data..., kwargs...))

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Convenience helpers
# -----------------------------------------------------------------------------

"""
    update_plankton_args(plankton_args::NamedTuple, group::Symbol; kwargs...) -> NamedTuple

Return a copy of `plankton_args` with the group entry `group` updated.

This is the primary user-facing helper for customizing plankton group arguments
(e.g. `:P`, `:Z`) without mutating the original container.

Behavior:
- Keyword arguments matching existing keys in `plankton_args.<group>` override
  those keys directly.
- Any remaining keyword arguments are interpreted as overrides for the group's
  `pft` (if present) and applied via `patch(::PFTSpecification; ...)`.
"""
function update_plankton_args(plankton_args::NamedTuple, group::Symbol; kwargs...)
    hasproperty(plankton_args, group) || throw(ArgumentError("update_plankton_args: group $(group) is not present in plankton_args"))
    spec = getproperty(plankton_args, group)
    spec isa NamedTuple || throw(ArgumentError("update_plankton_args: plankton_args.$(group) must be a NamedTuple, got $(typeof(spec))"))

    # Split kwargs into (1) direct overrides for the group spec and (2) pft overrides.
    direct_pairs = Pair{Symbol, Any}[]
    pft_pairs = Pair{Symbol, Any}[]
    for (k, v) in kwargs
        if hasproperty(spec, k)
            push!(direct_pairs, k => v)
        else
            push!(pft_pairs, k => v)
        end
    end

    patched_spec = isempty(direct_pairs) ? spec : patch(spec; direct_pairs...)

    if !isempty(pft_pairs)
        hasproperty(patched_spec, :pft) || throw(ArgumentError("update_plankton_args: plankton_args.$(group) has no `pft` field; cannot apply PFT overrides"))
        pft = getproperty(patched_spec, :pft)
        pft isa PFTSpecification || throw(ArgumentError("update_plankton_args: plankton_args.$(group).pft must be a PFTSpecification, got $(typeof(pft))"))
        patched_pft = patch(pft; pft_pairs...)
        patched_spec = merge(patched_spec, (; pft = patched_pft))
    end

    return merge(plankton_args, (; (group => patched_spec),))
end

"""
    update_biogeochem_args(args::BiogeochemistrySpecification; kwargs...) -> BiogeochemistrySpecification

Return a copy of `args` with overridden fields.

This is a thin, intent-revealing wrapper around `patch(args; ...)`.
"""
update_biogeochem_args(args::BiogeochemistrySpecification; kwargs...) = patch(args; kwargs...)

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

# Backwards-compatibility: keep the old name but route to the new API.
"""
    update_group(plankton_args::NamedTuple, group::Symbol; kwargs...) -> NamedTuple

Deprecated alias for `update_plankton_args`.
"""
update_group(plankton_args::NamedTuple, group::Symbol; kwargs...) = update_plankton_args(plankton_args, group; kwargs...)
