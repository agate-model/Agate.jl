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
            elseif v isa AbstractVector
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
    patch_registry_groups(registry::ParamRegistry, factory; kwargs...) -> ParamRegistry
    patch_registry_groups(registry::ParamRegistry, factory, patches::NamedTuple) -> ParamRegistry

Patch *group-level vector parameters* (stored as `GroupVec`) by updating only the
specified groups.

Each keyword in `kwargs...` must be the name of a **group-level vector** parameter in
`registry`. The corresponding value must be a `NamedTuple` mapping one or more group
symbols to scalar items (numbers, `Bool`, or `AbstractParamDef`).

Guarantees
----------
- Only the specified groups are updated; all other groups keep their existing values.
- Parameter names are validated strictly (typos throw).
- Group keys are validated strictly against the parameter's declared group set.
- `Dict` inputs are intentionally unsupported; use a `NamedTuple` like `(Z=...,)`.

Notes
-----
This function operates only on CPU-side registry metadata. It does **not** construct
any dense per-PFT vectors during patching.
"""
function patch_registry_groups(registry::ParamRegistry, factory, patches::NamedTuple)
    isempty(patches) && return registry

    new_specs = copy(registry.specs)
    for (param, patch) in pairs(patches)
        i = _lookup_index(registry, param)
        i == 0 && throw(ArgumentError("patch_registry_groups: parameter :$param is not present in this registry"))
        s = new_specs[i]

        s.shape === :vector || throw(ArgumentError(
            "patch_registry_groups: parameter :$param is not vector-shaped (shape=$(s.shape)). Only group-level vector parameters can be patched.",
        ))
        s.provider isa GroupVec || throw(ArgumentError(
            "patch_registry_groups: parameter :$param is not a group-level vector (provider=$(typeof(s.provider))). Use update_registry for strict replacement.",
        ))

        base = s.provider

        if patch isa AbstractDict
            throw(ArgumentError(
                "patch_registry_groups: group patches must be provided as a NamedTuple like (Z=...,). Dict inputs are intentionally unsupported.",
            ))
        end
        patch isa NamedTuple || throw(ArgumentError(
            "patch_registry_groups: patch for :$param must be a NamedTuple like (Z=...,). Got $(typeof(patch)).",
        ))

        # Validate patch keys (strict) and normalize values.
        patch_keys = Tuple(keys(patch))
        for g in patch_keys
            (g in base.groups) || throw(ArgumentError(
                "patch_registry_groups: unknown group :$g for parameter :$param. Expected a subset of groups=$(base.groups).",
            ))
        end

        # Normalize only the provided patch items; keep base items for others.
        N = length(base.groups)
        new_items = ntuple(i -> begin
            g = base.groups[i]
            if hasproperty(patch, g)
                _normalize_required_scalar_item(getproperty(patch, g), s.value_kind)
            else
                base.items[i]
            end
        end, N)

        gv = GroupVec{N}(base.groups, new_items)
        gv_norm = normalize_provider(:vector, gv, s.value_kind)
        new_specs[i] = _with_provider(s, gv_norm)
    end

    return ParamRegistry(new_specs, copy(registry.index))
end

function patch_registry_groups(registry::ParamRegistry, factory; kwargs...)
    isempty(kwargs) && return registry
    return patch_registry_groups(registry, factory, (; kwargs...))
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
