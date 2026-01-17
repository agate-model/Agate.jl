# This file is included into the `Agate.Parameters` module.

# ----------------------------------------------------------------------------
# Registry updates
# ----------------------------------------------------------------------------

"""
    update_registry(registry::ParamRegistry; kwargs...) -> ParamRegistry
    update_registry(registry::ParamRegistry, overrides::NamedTuple) -> ParamRegistry

Return a copy of `registry` with updated parameter providers.

Keys are validated to exist in the registry (strict-by-default) to catch typos early.
"""
@inline _with_provider(s::ParamSpec, provider) =
    ParamSpec(s.name, s.shape, s.value_kind, s.doc, provider)


@inline function _validate_complete_group_mapping(op::AbstractString, param::Symbol, groups, nt::NamedTuple)
    got = Tuple(keys(nt))
    missing = Symbol[]
    extra = Symbol[]
    for g in groups
        (g in got) || push!(missing, g)
    end
    for k in got
        (k in groups) || push!(extra, k)
    end
    if !isempty(missing) || !isempty(extra)
        msg = "$(op): group-level vector :$(param) requires a complete group mapping for groups=$(groups)."
        !isempty(missing) && (msg *= " Missing: $(missing).")
        !isempty(extra) && (msg *= " Unknown: $(extra).")
        msg *= " Provide a NamedTuple with exactly these keys."
        throw(ArgumentError(msg))
    end
    return nothing
end

function update_registry(registry::ParamRegistry, overrides::NamedTuple)
    isempty(overrides) && return registry

    new_specs = copy(registry.specs)
    for (k, v) in pairs(overrides)
        i = _lookup_index(registry, k)
        i == 0 && throw(ArgumentError("update_registry: parameter :$k is not present in this registry"))
        s = new_specs[i]

        new_provider = if s.shape === :vector && s.provider isa GroupVec
            # Group-level vectors are strict replacements.
            base = s.provider

            v === nothing && throw(ArgumentError(
                "update_registry: cannot unset group-level vector :$k. Provide a complete replacement for groups=$(base.groups).",
            ))

            gv = if v isa GroupVec
                v

            elseif v isa NamedTuple
                _validate_complete_group_mapping("update_registry", k, base.groups, v)
                GroupVec(base.groups, v)

            elseif v isa AbstractDict
                throw(ArgumentError(
                    "update_registry: :$k is a group-level vector parameter; provide a complete group NamedTuple matching groups=$(base.groups) (e.g. (P=1e-6, Z=1e-6)), or a GroupVec. Dict inputs are intentionally unsupported.",
                ))

            elseif v isa AbstractVector
                throw(ArgumentError(
                    "update_registry: :$k is a group-level vector parameter; provide a complete group NamedTuple matching groups=$(base.groups) (e.g. (P=1e-6, Z=1e-6)), or a GroupVec. A dense per-PFT vector is ambiguous for group-level parameters.",
                ))

            else
                throw(ArgumentError(
                    "update_registry: :$k is a group-level vector parameter; scalar broadcast is not supported. " *
                    "Provide a complete group NamedTuple matching groups=$(base.groups) (e.g. (P=1e-6, Z=1e-6)), or a GroupVec.",
                ))
            end

            gv_norm = normalize_provider(:vector, gv, s.value_kind)
            gv_norm.groups == base.groups || throw(ArgumentError(
                "update_registry: group-level vector :$k must use groups=$(base.groups), got groups=$(gv_norm.groups).",
            ))
            gv_norm

        elseif v === nothing
            nothing

        else
            normalize_provider(s.shape, v, s.value_kind)
        end

        new_specs[i] = _with_provider(s, new_provider)
    end

    # Providers only; parameter names do not change.
    return ParamRegistry(new_specs, copy(registry.index))
end

function update_registry(registry::ParamRegistry; kwargs...)
    isempty(kwargs) && return registry
    return update_registry(registry, (; kwargs...))
end
