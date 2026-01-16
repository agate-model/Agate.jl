# This file is included into the `Agate.Parameters` module.

# ----------------------------------------------------------------------------
# Registry updates
# ----------------------------------------------------------------------------

"""\
    update_registry(registry::ParamRegistry; kwargs...) -> ParamRegistry
    update_registry(registry::ParamRegistry, overrides::NamedTuple) -> ParamRegistry

Return a copy of `registry` with updated parameter providers.

Keys are validated to exist in the registry (strict-by-default) to catch typos early.
Use `extend_registry` to add new parameters explicitly.
"""
@inline _with_provider(s::ParamSpec, provider) =
    ParamSpec(s.name, s.shape, s.missing_policy, s.value_kind, s.doc, provider)

@inline function _merge_groupmap(a::VectorGroupMap, b::VectorGroupMap)
    keys = copy(a.keys)
    items = copy(a.items)
    for j in eachindex(b.keys)
        k = b.keys[j]
        it = b.items[j]

        idx = findfirst(==(k), keys)
        if idx === nothing
            push!(keys, k)
            push!(items, it)
        else
            items[idx] = it
        end
    end
    return VectorGroupMap(keys, items)
end

@inline function _patch_vector_provider(base, patch::VectorGroupMap)
    if base isa VectorGroupPatch
        # Compose patches: later updates override earlier ones for the same group.
        merged = _merge_groupmap(base.patch, patch)
        return VectorGroupPatch(base.base, merged)
    else
        return VectorGroupPatch(base, patch)
    end
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
                GroupVec(base.groups, v)
            elseif v isa AbstractDict
                throw(ArgumentError(
                    "update_registry: :$k is a group-level vector parameter; provide a scalar broadcast or a complete group NamedTuple matching groups=$(base.groups). Dict inputs are intentionally unsupported.",
                ))
            elseif v isa AbstractVector || v isa VectorGroupMap || v isa VectorGroupPatch
                throw(ArgumentError(
                    "update_registry: :$k is a group-level vector parameter; provide a scalar broadcast or a complete group NamedTuple matching groups=$(base.groups).",
                ))
            else
                # Scalar/Bool/allometric values broadcast across all groups.
                GroupVec{length(base.groups)}(base.groups, ntuple(_ -> v, length(base.groups)))
            end

            gv_norm = normalize_provider(:vector, gv, s.value_kind)
            gv_norm.groups == base.groups || throw(ArgumentError(
                "update_registry: group-level vector :$k must use groups=$(base.groups), got groups=$(gv_norm.groups).",
            ))
            gv_norm

        elseif v === nothing
            nothing

        elseif s.shape === :vector && (v isa NamedTuple || v isa VectorGroupMap)
            # Legacy per-group overrides are patches: update only the provided groups,
            # and delegate all other groups to the previous provider.
            patch = v isa VectorGroupMap ? v : normalize_provider(:vector, v, s.value_kind)
            _patch_vector_provider(s.provider, patch)

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

"""\
    extend_registry(registry::ParamRegistry, specs::ParamSpec...) -> ParamRegistry

Return a copy of `registry` with additional parameter specifications appended.

To avoid silent mistakes, extending with a key that already exists throws.
"""
function extend_registry(registry::ParamRegistry, specs::ParamSpec...)
    isempty(specs) && return registry

    new_specs = copy(registry.specs)
    new_index = copy(registry.index)
    for s in specs
        haskey(new_index, s.name) && throw(ArgumentError("extend_registry: parameter :$(s.name) already exists in registry"))
        push!(new_specs, s)
        new_index[s.name] = length(new_specs)
    end

    return ParamRegistry(new_specs, new_index)
end

"""\
    extend_registry(registry::ParamRegistry, specs::AbstractVector{ParamSpec}) -> ParamRegistry

Vector form of `extend_registry`.
"""
extend_registry(registry::ParamRegistry, specs::AbstractVector{ParamSpec}) =
    extend_registry(registry, specs...)
