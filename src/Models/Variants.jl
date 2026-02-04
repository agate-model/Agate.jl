"""Model-variant scaffolding.

As models grow, we want to avoid copying entire modules just to represent a
paper-specific configuration or an internal experiment.

This file provides a small registry for **variants**: lightweight construction-time
recipes that bundle together a factory, default community structure, default
dynamics, and optional parameter/interaction overrides.

Naming
------
Variants are identified by a `ModelId(family, citation, tag)` where:
- `family`: the model family (e.g. `:DARWIN`, `:NiPiZD`).
- `citation`: a citation key (e.g. `:citation2026`).
- `tag`: a stable label for a concrete recipe within that citation (e.g. `:A`, `:B`,
  `:submission`, `:accepted`). We intentionally avoid implying chronology; use any
  tag scheme that remains meaningful for your workflow.

This scaffolding is intentionally small. It is *not* a full experiment-management
system; it only helps keep model-family code clean as variants accumulate.
"""

using ..Utils: AbstractBGCFactory
using ..Constructor: construct_factory

"""Stable identifier for a model variant.

Examples
--------
```julia
ModelId(:DARWIN, :citation2026, :A)
ModelId(:DARWIN, :citation2026, :submission)
```
"""
struct ModelId
    family::Symbol
    citation::Symbol
    tag::Symbol
end

Base.string(id::ModelId) = string(id.family, "/", id.citation, "/", id.tag)

"""Construction-time container for a model variant.

`VariantSpec` stores the *inputs* to `Agate.Constructor.construct_factory`.
"""
struct VariantSpec{
    F<:AbstractBGCFactory,
    PD<:NamedTuple,
    BD<:NamedTuple,
    C<:NamedTuple,
    R<:NamedTuple,
    PG<:NamedTuple,
    AF<:Tuple,
    P<:NamedTuple,
    I,
}
    id::ModelId
    factory::F
    plankton_dynamics::PD
    biogeochem_dynamics::BD
    community::C
    roles::R
    parameter_groups::PG
    auxiliary_fields::AF
    parameters::P
    interaction_overrides::I # `Nothing` or a `NamedTuple`
end

"""Internal registry of known variants.

Values are builder functions that accept keyword arguments and return a `VariantSpec`.
"""
const _VARIANT_REGISTRY = Dict{ModelId,Function}()

"""Register a variant builder.

The builder must return a `VariantSpec`.
"""
function register_variant(id::ModelId, builder::Function)
    _VARIANT_REGISTRY[id] = builder
    return id
end

"""Return the list of registered variants.

You may filter by `family` and/or `citation`.
"""
function list_variants(;
    family::Union{Nothing,Symbol}=nothing, citation::Union{Nothing,Symbol}=nothing
)
    ids = collect(keys(_VARIANT_REGISTRY))
    if !isnothing(family)
        ids = filter(id -> id.family == family, ids)
    end
    if !isnothing(citation)
        ids = filter(id -> id.citation == citation, ids)
    end
    sort!(ids; by=string)
    return ids
end

"""Construct a `VariantSpec` from the registry."""
function variant(id::ModelId; kwargs...)
    builder = get(_VARIANT_REGISTRY, id, nothing)
    isnothing(builder) && throw(
        ArgumentError(
            "unknown variant $(string(id)); available variants: $(join(string.(list_variants()), ", "))",
        ),
    )
    return builder(; kwargs...)
end

"""Convenience overload: `variant(:DARWIN, :citation2026, :A; kwargs...)`."""
variant(family::Symbol, citation::Symbol, tag::Symbol; kwargs...) =
    variant(ModelId(family, citation, tag); kwargs...)

"""Construct a model instance from a `VariantSpec`.

This is a thin convenience wrapper around `Agate.Constructor.construct_factory`.
"""
function construct(
    spec::VariantSpec;
    parameters::NamedTuple=(;),
    interaction_overrides::Union{Nothing,NamedTuple}=nothing,
    roles=nothing,
    parameter_groups=nothing,
    auxiliary_fields=nothing,
    kwargs...,
)
    # Merge runtime overrides on top of the variant defaults.
    params = merge(spec.parameters, parameters)

    inter = spec.interaction_overrides
    if isnothing(inter)
        inter = interaction_overrides
    elseif isnothing(interaction_overrides)
        # Keep the spec's interaction overrides.
    else
        inter = merge(inter, interaction_overrides)
    end

    return construct_factory(
        spec.factory;
        plankton_dynamics=spec.plankton_dynamics,
        biogeochem_dynamics=spec.biogeochem_dynamics,
        community=spec.community,
        parameters=params,
        interaction_overrides=inter,
        roles=isnothing(roles) ? spec.roles : roles,
        parameter_groups=isnothing(parameter_groups) ? spec.parameter_groups : parameter_groups,
        auxiliary_fields=isnothing(auxiliary_fields) ? spec.auxiliary_fields : auxiliary_fields,
        kwargs...,
    )
end
