using OceanBioME: BoxModelGrid

import ...Configuration
import ...Construction
import ...Factories

export construct

"""
    construct(; kw...) -> bgc

Construct a size-structured NiPiZD ecosystem model.

The NiPiZD model contains two plankton groups: phytoplankton (`P`) and zooplankton (`Z`),
each represented by `n_phyto` and `n_zoo` size classes.

In addition to plankton, the default NiPiZD factory includes idealized nutrient (`N`) and
detritus (`D`) cycling. The returned biogeochemistry instance includes a photosynthetically
active radiation (PAR) auxiliary field.

During construction, plankton size (diameter) is used to resolve trait-based parameter
vectors and interaction matrices (e.g. palatability and assimilation efficiency). You
may override interaction matrices explicitly with `palatability_matrix` and/or
`assimilation_matrix`.

Keywords
--------
- `n_phyto=2`, `n_zoo=2`: number of phytoplankton and zooplankton size classes
- `phyto_diameters=(2, 10, :log_splitting)`: diameter specification for phytoplankton
- `zoo_diameters=(20, 100, :linear_splitting)`: diameter specification for zooplankton
- `parameters=(;)`: parameter overrides (validated against the NiPiZD parameter set)
- `palatability_matrix=nothing`, `assimilation_matrix=nothing`: optional interaction matrix overrides
- `grid=BoxModelGrid()`: grid used for precision/architecture inference and sinking velocity fields
- `arch=nothing`: override the architecture (usually inferred from `grid`)
- `sinking_tracers=nothing`: sinking speed overrides, e.g. `(D = 2/day, P1 = 0.1/day, ...)`
- `open_bottom=true`: whether sinking tracers leave the domain

Returns
-------
An `Oceananigans.Biogeochemistry.AbstractContinuousFormBiogeochemistry` instance.

Example
-------
```julia
using Agate.Models: NiPiZD

bgc = NiPiZD.construct()
```
"""
function construct(;
    n_phyto::Int=2,
    n_zoo::Int=2,
    phyto_diameters=(2, 10, :log_splitting),
    zoo_diameters=(20, 100, :linear_splitting),
    parameters::NamedTuple=(;),
    palatability_matrix=nothing,
    assimilation_matrix=nothing,
    grid=BoxModelGrid(),
    arch=nothing,
    sinking_tracers=nothing,
    open_bottom::Bool=true,
)
    n_phyto >= 1 || throw(ArgumentError("n_phyto must be >= 1"))
    n_zoo >= 1 || throw(ArgumentError("n_zoo must be >= 1"))

    factory = NiPiZDFactory()

    base = Factories.default_community(factory)
    community = Configuration.build_plankton_community(
        base;
        n=(Z=n_zoo, P=n_phyto),
        diameters=(
            Z=Configuration.diameter_specification(zoo_diameters),
            P=Configuration.diameter_specification(phyto_diameters),
        ),
    )

    # Interaction overrides (optional).
    #
    # We forward overrides through the model-agnostic constructor as a `NamedTuple`.
    # Interaction overrides are data-only: values must be explicit matrices.
    pairs = Pair{Symbol,Any}[]
    palatability_matrix !== nothing &&
        push!(pairs, :palatability_matrix => palatability_matrix)
    assimilation_matrix !== nothing &&
        push!(pairs, :assimilation_matrix => assimilation_matrix)

    interaction_overrides = isempty(pairs) ? nothing : (; pairs...)

    interaction_roles = (consumers=(:Z,), prey=(:P,))

    return Construction.construct_factory(
        factory;
        community=community,
        parameters=parameters,
        interaction_roles=interaction_roles,
        auxiliary_fields=(:PAR,),
        interaction_overrides=interaction_overrides,
        arch=arch,
        sinking_tracers=sinking_tracers,
        grid=grid,
        open_bottom=open_bottom,
    )
end
