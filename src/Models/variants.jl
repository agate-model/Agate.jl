using ..Factories: AbstractBGCFactory
using ..Construction: construct_factory_with_context
using ..Traceability: finalize_construction_manifest

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

`VariantSpec` stores the *inputs* to `Agate.Construction.construct_factory`.
"""
struct VariantSpec{
    F<:AbstractBGCFactory,
    PD<:NamedTuple,
    BD<:NamedTuple,
    C<:NamedTuple,
    R<:NamedTuple,
    AF<:Tuple,
    P<:NamedTuple,
    I,
}
    id::ModelId
    factory::F
    plankton_dynamics::PD
    biogeochem_dynamics::BD
    community::C
    interaction_roles::R
    auxiliary_fields::AF
    parameters::P
    interaction_overrides::I # `Nothing` or a `NamedTuple`
end

"""Internal registry of known variants.

Values are builder functions that accept keyword arguments and return a `VariantSpec`.
"""
const VARIANT_REGISTRY = Dict{ModelId,Function}()

"""Register a variant builder.

The builder must return a `VariantSpec`.
"""
function register_variant(id::ModelId, builder::Function)
    VARIANT_REGISTRY[id] = builder
    return id
end

"""Return the list of registered variants.

You may filter by `family` and/or `citation`.
"""
function list_variants(;
    family::Union{Nothing,Symbol}=nothing, citation::Union{Nothing,Symbol}=nothing
)
    ids = collect(keys(VARIANT_REGISTRY))
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
    builder = get(VARIANT_REGISTRY, id, nothing)
    isnothing(builder) && throw(
        ArgumentError(
            "unknown variant $(string(id)); available variants: $(join(string.(list_variants()), ", "))",
        ),
    )
    return builder(; kwargs...)
end

"""Convenience overload: `variant(:DARWIN, :citation2026, :A; kwargs...)`."""
function variant(family::Symbol, citation::Symbol, tag::Symbol; kwargs...)
    variant(ModelId(family, citation, tag); kwargs...)
end

"""Construct a model instance from a `VariantSpec`.

This is a thin convenience wrapper around `Agate.Construction.construct_factory`.
"""
function construct(
    spec::VariantSpec;
    parameters::NamedTuple=(;),
    interaction_overrides::Union{Nothing,NamedTuple}=nothing,
    auxiliary_fields=nothing,
    kwargs...,
)
    bgc, _ = _construct_variant_with_context(
        spec;
        parameters=parameters,
        interaction_overrides=interaction_overrides,
        auxiliary_fields=auxiliary_fields,
        kwargs...,
    )
    return bgc
end

"""Construct a model instance and return `(bgc, manifest)`.

The returned `bgc` is the same runtime object returned by `construct(spec; ...)`.
The manifest is generated during construction from the same resolved state and is
kept separate from `AgateBGC` so the runtime object stays GPU-friendly.
"""
function construct_with_manifest(
    spec::VariantSpec;
    parameters::NamedTuple=(;),
    interaction_overrides::Union{Nothing,NamedTuple}=nothing,
    auxiliary_fields=nothing,
    kwargs...,
)
    bgc, context = _construct_variant_with_context(
        spec;
        parameters=parameters,
        interaction_overrides=interaction_overrides,
        auxiliary_fields=auxiliary_fields,
        kwargs...,
    )
    return bgc, finalize_construction_manifest(context)
end

function _construct_variant_with_context(
    spec::VariantSpec;
    parameters::NamedTuple=(;),
    interaction_overrides::Union{Nothing,NamedTuple}=nothing,
    auxiliary_fields=nothing,
    kwargs...,
)
    params = merge(spec.parameters, parameters)

    inter = spec.interaction_overrides
    if isnothing(inter)
        inter = interaction_overrides
    elseif isnothing(interaction_overrides)
        # Keep the spec's interaction overrides.
    else
        inter = merge(inter, interaction_overrides)
    end

    effective_auxiliary_fields = if isnothing(auxiliary_fields)
        spec.auxiliary_fields
    else
        auxiliary_fields
    end

    bgc, resolved = construct_factory_with_context(
        spec.factory;
        plankton_dynamics=spec.plankton_dynamics,
        biogeochem_dynamics=spec.biogeochem_dynamics,
        community=spec.community,
        parameters=params,
        interaction_overrides=inter,
        interaction_roles=spec.interaction_roles,
        auxiliary_fields=effective_auxiliary_fields,
        kwargs...,
    )

    context = (
        model=(
            id=string(spec.id),
            family=string(spec.id.family),
            citation=string(spec.id.citation),
            tag=string(spec.id.tag),
        ),
        resolved=resolved,
    )

    return bgc, context
end
