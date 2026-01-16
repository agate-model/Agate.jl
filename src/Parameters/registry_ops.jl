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
    # Merge two group maps, with `b` overriding entries in `a` for duplicate keys.
    keys = copy(a.keys)
    items = copy(a.items)
    for (k, item) in zip(b.keys, b.items)
        j = findfirst(==(k), keys)
        if j === nothing
            push!(keys, k)
            push!(items, item)
        else
            items[j] = item
        end
    end
    return VectorGroupMap(keys, items)
end

function update_registry(registry::ParamRegistry, overrides::NamedTuple)
    isempty(overrides) && return registry

    new_specs = copy(registry.specs)
    for (k, v) in pairs(overrides)
        i = _lookup_index(registry, k)
        i == 0 && throw(ArgumentError("update_registry: parameter :$k is not present in this registry"))
        s = new_specs[i]
        new_provider = if v === nothing
            nothing
        elseif s.shape === :vector && v isa NamedTuple
            patch = normalize_provider(:vector, v, s.value_kind)
            patch isa VectorGroupMap || throw(ArgumentError("Internal error: vector group map normalization failed for :$k"))
            base = s.provider
            if base isa VectorGroupPatch
                merged = _merge_groupmap(base.patch, patch)
                VectorGroupPatch(base.base, merged)
            else
                VectorGroupPatch(base, patch)
            end
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