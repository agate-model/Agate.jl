"""Agate.Constructor model specifications.

`ModelSpec` is an internal construction-time container used by Agate's public model
constructors (for example `NiPiZD.construct` and `DARWIN.construct`).

It exists solely to reduce boilerplate in model-specific wrappers by bundling a
factory with sensible defaults (community + dynamics), while leaving overrides as
plain `NamedTuple`s.
"""

using ..Utils: AbstractBGCFactory
using ..Interface:
    default_plankton_dynamics, default_biogeochem_dynamics, default_community

"""A lightweight wrapper around an `AbstractBGCFactory`."""
struct ModelSpec{F<:AbstractBGCFactory}
    factory::F
end

"""Build a `(Z,P)` community from a base community spec.

This is used by both NiPiZD and DARWIN, which share the same plankton group symbols.
"""
function build_ZP_community(
    base::NamedTuple; n_zoo::Int, n_phyto::Int, zoo_diameters, phyto_diameters
)
    Z = (; base.Z..., n=n_zoo, diameters=zoo_diameters)
    P = (; base.P..., n=n_phyto, diameters=phyto_diameters)
    return (Z=Z, P=P)
end

"""Construct a model instance from `spec`.

Keywords
--------
- `plankton_dynamics`: `NamedTuple` mapping group symbols (e.g. `:Z`, `:P`) to builder functions.
- `biogeochem_dynamics`: `NamedTuple` mapping non-plankton tracer symbols to builder functions.
- `community`: `NamedTuple` describing plankton size structure.
- `parameters`: `NamedTuple` of parameter overrides.
- `interaction_overrides`: optional interaction-related overrides (typically matrices). Values may be matrices or provider functions.
- `roles`: optional `NamedTuple` mapping group symbols to role symbols.
- `arch`: architecture specification.
- `sinking_tracers`: collection of tracer symbols that sink and their sinking velocities.
- `grid`: grid specification.
- `open_bottom`: boolean indicating if the bottom boundary is open.

All other keywords are forwarded to the model-agnostic constructor.
"""
function construct_factory(
    spec::ModelSpec;
    plankton_dynamics=default_plankton_dynamics(spec.factory),
    biogeochem_dynamics=default_biogeochem_dynamics(spec.factory),
    community=default_community(spec.factory),
    parameters::NamedTuple=(;),
    interaction_overrides::Union{Nothing,NamedTuple}=nothing,
    roles::Union{Nothing,NamedTuple}=nothing,
    arch=nothing,
    sinking_tracers=nothing,
    grid=nothing,
    open_bottom::Bool=true,
)
    return construct_factory(
        spec.factory;
        plankton_dynamics=plankton_dynamics,
        biogeochem_dynamics=biogeochem_dynamics,
        community=community,
        parameters=parameters,
        interaction_overrides=interaction_overrides,
        roles=roles,
        arch=arch,
        sinking_tracers=sinking_tracers,
        grid=grid,
        open_bottom=open_bottom,
    )
end
