"""Public constructor for the NiPiZD ecosystem model.

`NiPiZD.construct` builds an Oceananigans/OceanBioME-compatible biogeochemistry
instance for a size-structured plankton model with:

- two plankton groups: phytoplankton (`P`) and zooplankton (`Z`)
- two non-plankton tracers: dissolved inorganic nutrient (`N`) and detritus (`D`)

The constructor exposes a small set of keywords for community structure,
dynamics, parameter overrides, and interaction-matrix overrides.
"""

using OceanBioME: BoxModelGrid

using .Tracers:
    nutrient_default,
    detritus_default,
    phytoplankton_default,
    zooplankton_default

import ...Utils
import ...Constructor
import ...FactoryInterface

export construct

@inline function _call_matrix_provider(f::Function, ctx)
    try
        return f(ctx.diameters, ctx.group_symbols)
    catch err
        err isa MethodError || rethrow()
        return f(ctx)
    end
end

@inline _materialize_matrix_override(x, ctx) = x isa Function ? _call_matrix_provider(x, ctx) : x

function _nipizd_interactions(palatability_matrix, assimilation_matrix)
    palatability_matrix === nothing && assimilation_matrix === nothing && return nothing

    has_provider = (palatability_matrix isa Function) || (assimilation_matrix isa Function)

    if has_provider
        return function (ctx)
            overrides = (;)

            if palatability_matrix !== nothing
                overrides = merge(overrides, (; palatability_matrix = _materialize_matrix_override(palatability_matrix, ctx)))
            end

            if assimilation_matrix !== nothing
                overrides = merge(overrides, (; assimilation_matrix = _materialize_matrix_override(assimilation_matrix, ctx)))
            end

            return overrides
        end
    end

    overrides = (;)
    palatability_matrix !== nothing && (overrides = merge(overrides, (; palatability_matrix = palatability_matrix)))
    assimilation_matrix !== nothing && (overrides = merge(overrides, (; assimilation_matrix = assimilation_matrix)))

    return overrides
end

"""
    construct(; kw...) -> bgc

Construct a NiPiZD biogeochemistry instance.

Keywords
--------
- `n_phyto=2`, `n_zoo=2`: number of phytoplankton and zooplankton size classes
- `phyto_diameters=(2, 10, :log_splitting)`: diameter specification for phytoplankton
- `zoo_diameters=(20, 100, :linear_splitting)`: diameter specification for zooplankton
- `nutrient_dynamics`, `detritus_dynamics`, `phyto_dynamics`, `zoo_dynamics`: dynamics builders
- `parameters=(;)`: parameter overrides
- `palatability_matrix=nothing`: matrix (or provider) overriding `:palatability_matrix`
- `assimilation_matrix=nothing`: matrix (or provider) overriding `:assimilation_matrix`
- `grid=BoxModelGrid()`: grid used for precision/architecture inference and sinking velocity fields
- `arch=nothing`: override the architecture (usually inferred from `grid`)
- `sinking_tracers=nothing`: sinking speed overrides, e.g. `(D = 2/ day, P1 = 0.1/day, ...)`
- `open_bottom=true`: whether sinking tracers leave the domain

Matrix providers
----------------
If `palatability_matrix` or `assimilation_matrix` is a function, it is called during
construction as either:

- `f(diameters, group_symbols)`
- `f(ctx)` (fallback for single-argument providers)

Returns
-------
An `Oceananigans.Biogeochemistry.AbstractContinuousFormBiogeochemistry` instance.
"""
function construct(;
    n_phyto::Int = 2,
    n_zoo::Int = 2,
    phyto_diameters = (2, 10, :log_splitting),
    zoo_diameters = (20, 100, :linear_splitting),
    nutrient_dynamics = nutrient_default,
    detritus_dynamics = detritus_default,
    phyto_dynamics = phytoplankton_default,
    zoo_dynamics = zooplankton_default,
    parameters::NamedTuple = (;),
    palatability_matrix = nothing,
    assimilation_matrix = nothing,
    grid = BoxModelGrid(),
    arch = nothing,
    sinking_tracers = nothing,
    open_bottom::Bool = true,
)
    n_phyto >= 1 || throw(ArgumentError("n_phyto must be >= 1"))
    n_zoo >= 1 || throw(ArgumentError("n_zoo must be >= 1"))

    factory = NiPiZDFactory()
    spec = Constructor.ModelSpec(factory)

    base = FactoryInterface.default_community(factory)
    community = Constructor.build_ZP_community(
        base;
        n_zoo = n_zoo,
        n_phyto = n_phyto,
        zoo_diameters = Utils.diameter_specification(zoo_diameters),
        phyto_diameters = Utils.diameter_specification(phyto_diameters),
    )

    plankton_dynamics = (Z = zoo_dynamics, P = phyto_dynamics)
    biogeochem_dynamics = (N = nutrient_dynamics, D = detritus_dynamics)

    interactions = _nipizd_interactions(palatability_matrix, assimilation_matrix)

    return Constructor.construct(
        spec;
        plankton_dynamics = plankton_dynamics,
        biogeochem_dynamics = biogeochem_dynamics,
        community = community,
        parameters = parameters,
        interactions = interactions,
        arch = arch,
        sinking_tracers = sinking_tracers,
        grid = grid,
        open_bottom = open_bottom,
    )
end
