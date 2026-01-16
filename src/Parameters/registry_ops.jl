# This file is included into the `Agate.Parameters` module.

# ----------------------------------------------------------------------------
# Registry updates
# ----------------------------------------------------------------------------

"""\
    update_registry(registry::ParamRegistry; kwargs...) -> ParamRegistry

Return a copy of `registry` with updated parameter providers.

Keys are validated to exist in the registry (strict-by-default) to catch typos early.
Use `extend_registry` to add new parameters explicitly.
"""
function update_registry(registry::ParamRegistry; kwargs...)
    isempty(kwargs) && return registry

    overrides = (; kwargs...)
    new_specs = copy(registry.specs)

    for (k, v) in pairs(overrides)
        i = _lookup_index(registry, k)
        i == 0 && throw(ArgumentError("update_registry: parameter $(k) is not present in this registry"))
        s = new_specs[i]
        new_provider = v === nothing ? nothing : normalize_provider(s.shape, v)
        new_specs[i] = ParamSpec(s.name, s.shape, s.missing_policy, s.value_kind, s.doc, new_provider)
    end

    return ParamRegistry(new_specs)
end

"""\
    extend_registry(registry::ParamRegistry, specs::ParamSpec...) -> ParamRegistry

Return a copy of `registry` with additional parameter specifications appended.

To avoid silent mistakes, extending with a key that already exists throws.
"""
function extend_registry(registry::ParamRegistry, specs::ParamSpec...)
    isempty(specs) && return registry

    for s in specs
        haskey(registry.index, s.name) && throw(ArgumentError("extend_registry: parameter $(s.name) already exists in registry"))
    end

    new_specs = copy(registry.specs)
    append!(new_specs, specs)
    return ParamRegistry(new_specs)
end
