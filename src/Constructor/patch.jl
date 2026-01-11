

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
# Convenience helpers
# -----------------------------------------------------------------------------

"""
    update_group(plankton_args::NamedTuple, group::Symbol; kwargs...) -> NamedTuple

Convenience wrapper for updating a single plankton group entry in `plankton_args`.

Equivalent to:

```julia
patch(plankton_args; group = (kwargs...,))
```

but allows `group` to be specified programmatically.
"""
function update_group(plankton_args::NamedTuple, group::Symbol; kwargs...)
    hasproperty(plankton_args, group) || throw(ArgumentError("update_group: group $(group) is not present in plankton_args"))
    spec = getproperty(plankton_args, group)
    spec isa NamedTuple || throw(ArgumentError("update_group: plankton_args.$(group) must be a NamedTuple, got $(typeof(spec))"))

    patched_spec = patch(spec; kwargs...)
    return merge(plankton_args, (; (group => patched_spec),))
end
