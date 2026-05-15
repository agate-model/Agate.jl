using OceanBioME: BoxModelGrid

import ...Configuration
import ...Construction
import ...Factories

export construct

"""
    construct(; kw...) -> bgc

Construct a size-structured NiPiZD ecosystem model.

The NiPiZD model contains two plankton groups: phytoplankton (`P`) and zooplankton (`Z`),
with size classes defined by their size-structure inputs.

In addition to plankton, the default NiPiZD factory includes idealized nutrient (`N`) and
detritus (`D`) cycling. The returned biogeochemistry instance includes a photosynthetically
active radiation (PAR) auxiliary field.

During construction, plankton size (diameter) is used to resolve trait-based parameter
vectors and interaction matrices (e.g. palatability and assimilation efficiency). You
may override interaction matrices explicitly with `palatability_matrix` and/or
`assimilation_matrix`.

Size-structure inputs may be a NamedTuple range, for example
`(n=3, min_esd=1, max_esd=10, splitting=:log_splitting)`, or an explicit
diameter vector such as `[1.0, 3.2, 10.0]`.

Keywords
--------
- `phyto_size_structure=(n=2, min_esd=2, max_esd=10, splitting=:log_splitting)`: phytoplankton size structure
- `zoo_size_structure=(n=2, min_esd=20, max_esd=100, splitting=:linear_splitting)`: zooplankton size structure
- `parameters=(;)`: parameter overrides (validated against the NiPiZD parameter set)
- `palatability_matrix=nothing`: optional palatability matrix override. Must be an explicit rectangular matrix sized to the canonical interaction axes `(n_consumer, n_prey)`.
- `assimilation_matrix=nothing`: optional assimilation matrix override. Must be an explicit rectangular matrix sized to the canonical interaction axes `(n_consumer, n_prey)`.
- `grid=BoxModelGrid()`: grid used for architecture inference and default scalar-type selection
- `scalar_type=nothing`: explicit runtime scalar type. When omitted, construction uses `eltype(grid)` or `Float64` if no grid is supplied
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
    phyto_size_structure=(n=2, min_esd=2, max_esd=10, splitting=:log_splitting),
    zoo_size_structure=(n=2, min_esd=20, max_esd=100, splitting=:linear_splitting),
    parameters::NamedTuple=(;),
    palatability_matrix=nothing,
    assimilation_matrix=nothing,
    grid=BoxModelGrid(),
    scalar_type=nothing,
    arch=nothing,
    sinking_tracers=nothing,
    open_bottom::Bool=true,
)
    factory = NiPiZDFactory()

    base = Factories.default_community(factory)
    community = Configuration.build_plankton_community(
        base; diameters=(Z=zoo_size_structure, P=phyto_size_structure)
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
        scalar_type=scalar_type,
        open_bottom=open_bottom,
    )
end
