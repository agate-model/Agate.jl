"""Agate.Constructor model specifications.

`ModelSpec` is an internal construction-time container used by Agate's public model
constructors (for example `NiPiZD.construct` and `DARWIN.construct`).

It implements a small, explicit **spec + merge** workflow:

1. Start from model defaults (community, dynamics, parameter registry).
2. Merge user overrides into those defaults.
3. Delegate to the model-agnostic constructor (`construct(factory; ...)`).

`ModelSpec` is intentionally minimal: it reduces boilerplate for model-specific wrappers
without introducing an extensible public "factory API" surface.
"""

using ..Utils: AbstractBGCFactory
using ..FactoryInterface:
    default_plankton_dynamics,
    default_biogeochem_dynamics,
     default_community

import ..Parameters

"""A lightweight wrapper around an `AbstractBGCFactory`.

A `ModelSpec` is not itself a factory. It exists so public constructors can share a
single merge-based implementation.
"""
struct ModelSpec{F<:AbstractBGCFactory}
    factory::F
end

ModelSpec(factory::AbstractBGCFactory) = ModelSpec{typeof(factory)}(factory)

@inline _isempty(nt::NamedTuple) = isempty(keys(nt))

"""Build a `(Z,P)` community from a base community spec.

This is used by both NiPiZD and DARWIN, which share the same plankton group symbols.
"""
function build_ZP_community(base::NamedTuple; n_zoo::Int, n_phyto::Int, zoo_diameters, phyto_diameters)
    Z = (; base.Z..., n=n_zoo, diameters=zoo_diameters)
    P = (; base.P..., n=n_phyto, diameters=phyto_diameters)
    return (Z = Z, P = P)
end

"""Construct a model instance from `spec` by merging explicit overrides.

Keywords
--------
- `plankton_dynamics`: `NamedTuple` mapping group symbols (e.g. `:Z`, `:P`) to builder functions.
- `biogeochem_dynamics`: `NamedTuple` mapping non-plankton tracer symbols (e.g. `:N`, `:DIC`) to builder functions.
- `community`: `NamedTuple` describing plankton size structure.
- `parameters`: `NamedTuple` of parameter overrides (validated against the model's registry).
- `interactions`: optional interaction-related overrides (typically matrices).

All other keywords are forwarded to the model-agnostic constructor.
"""
function construct(
    spec::ModelSpec;
    plankton_dynamics = default_plankton_dynamics(spec.factory),
    biogeochem_dynamics = default_biogeochem_dynamics(spec.factory),
    community = default_community(spec.factory),
    parameters::NamedTuple = (;),
    registry = nothing,
    interactions::Union{Nothing,NamedTuple,Function} = nothing,
    arch = nothing,
    sinking_tracers = nothing,
    grid = nothing,
    open_bottom::Bool = true,
)
    reg = isnothing(registry) ? Parameters.parameter_registry(spec.factory) : registry

    if !_isempty(parameters)
        reg = Parameters.update_registry(reg; parameters...)
    end

    return construct(
        spec.factory;
        plankton_dynamics = plankton_dynamics,
        biogeochem_dynamics = biogeochem_dynamics,
        community = community,
        registry = reg,
        interactions = interactions,
        arch = arch,
        sinking_tracers = sinking_tracers,
        grid = grid,
        open_bottom = open_bottom,
    )
end
