using OceanBioME: BoxModelGrid

import ...Configuration
import ...Construction
import ...Factories

using ..: default_model_manifest_context, default_model_recipe, finalize_construction_manifest

export construct, construct_with_manifest

function _construction_inputs(;
    phyto_size_structure=(n=2, min_esd=1.5, max_esd=20.0, splitting=:log_splitting),
    zoo_size_structure=(n=2, min_esd=20.0, max_esd=100.0, splitting=:log_splitting),
    parameters::NamedTuple=(;),
    scalar_type=nothing,
    palatability_matrix=nothing,
    assimilation_matrix=nothing,
    grid=BoxModelGrid(),
    arch=nothing,
    sinking_tracers=nothing,
    open_bottom::Bool=true,
)
    factory = DarwinFactory()

    base = Factories.default_community(factory)
    community = Configuration.build_plankton_community(
        base; diameters=(Z=zoo_size_structure, P=phyto_size_structure)
    )

    pairs = Pair{Symbol,Any}[]
    palatability_matrix !== nothing &&
        push!(pairs, :palatability_matrix => palatability_matrix)
    assimilation_matrix !== nothing &&
        push!(pairs, :assimilation_matrix => assimilation_matrix)

    interaction_overrides = isempty(pairs) ? nothing : (; pairs...)
    interaction_roles = (consumers=(:Z,), prey=(:P,))

    return (
        factory=factory,
        kwargs=(;
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
        ),
    )
end

"""
    construct(; kw...) -> bgc
    construct_with_manifest(; kw...) -> bgc, manifest

Construct a simplified DARWIN-like, size-structured ecosystem model.

!!! info
    This model is in active development and has not been validated against `MITgcm-DARWIN`.

!!! formulation
    TRACERS:

    ∂t cⱼ = ``Uⱼ``DIC - ``Mⱼ`` + ``Gⱼ`` - ``gⱼ``

    ∂t DIC = ∑(``Uⱼ`` DIC) + ``R``DOC + ``R``POC

    ∂t DIN = ∑(``Uⱼ``DIC * ``Qⱼ``N)  + ``R``DON + ``R``PON

    ∂t PO4 = ∑(``Uⱼ``DIC * ``Qⱼ``P)  + ``R``DOP + ``R``POP

    ∂t DOC = ∑(``Mⱼ``DOC) + ``g``DOC - ``R``DOC

    ∂t DON = ∑(``Mⱼ``DOC * ``Qⱼ``N) + ``g``DON - ``R``DON

    ∂t DOP = ∑(``Mⱼ``DOC * ``Qⱼ``P) + ``g``DOP - ``R``DOP

    ∂t POC = ∑(``Mⱼ``POC) + ``g``POC - ``R``POC

    ∂t PON = ∑(``Mⱼ``POC * ``Qⱼ``N) + ``g``PON - ``R``PON

    ∂t POP = ∑(``Mⱼ``POC * ``Qⱼ``P) + ``g``POP - ``R``POP

    where:
    - ``U`` = uptake
    - ``R`` = remineralization
    - ``M`` = mortality
    - ``g, G`` = grazing losses and gains
    - ``Q`` = plankton elemental ratios

    TRAITS:

    μmax, KR, gmax = a*Volume^b

    palat = η/(1+(``ratio``-``opt``)^2)^σ

    where:
    - μmax = maximum photosynthetic growth
    - KR = nutrient half saturation
    - gmax = maximum predation rate
    - palat = palatability
    - ``ratio`` = predator to prey size ratio (diameter)
    - ``opt`` = predator to prey size optimum (diameter)
    - η = prey protection
    - σ = predator specificity

The DARWIN model contains two plankton groups: phytoplankton (`P`) and zooplankton (`Z`),
with size classes defined by their size-structure inputs.

In addition to plankton, the default DARWIN factory includes elemental cycling tracers:
dissolved inorganic carbon (DIC), dissolved inorganic nitrogen (DIN), phosphate (PO4),
dissolved organic matter (DOC, DON, DOP), and particulate organic matter (POC, PON, POP).

During construction, plankton size (diameter) is used to resolve trait-based parameter
vectors and interaction matrices (e.g. palatability and assimilation efficiency). You
may override interaction matrices explicitly with `palatability_matrix` and/or
`assimilation_matrix`.

Size-structure inputs may be a NamedTuple range, for example
`(n=3, min_esd=1, max_esd=10, splitting=:log_splitting)`, or an explicit
diameter vector such as `[1.0, 3.2, 10.0]`.

The returned biogeochemistry instance includes a photosynthetically active radiation (PAR)
auxiliary field.

Keywords
--------
- `phyto_size_structure=(n=2, min_esd=1.5, max_esd=20.0, splitting=:log_splitting)`: phytoplankton size structure
- `zoo_size_structure=(n=2, min_esd=20.0, max_esd=100.0, splitting=:log_splitting)`: zooplankton size structure
- `parameters=(;)`: parameter overrides (validated against the DARWIN parameter set). Vector parameters may be supplied positionally or as partial NamedTuple overrides keyed by plankton class name.
- `scalar_type=nothing`: explicit runtime scalar type. When omitted, construction uses `eltype(grid)` or `Float64` if no grid is supplied
- `palatability_matrix=nothing`: optional palatability matrix override. Must be an explicit rectangular matrix sized to the canonical interaction axes `(n_consumer, n_prey)`.
- `assimilation_matrix=nothing`: optional assimilation matrix override. Must be an explicit rectangular matrix sized to the canonical interaction axes `(n_consumer, n_prey)`

- `grid=BoxModelGrid()`: grid used for precision/architecture inference and sinking velocity fields
- `arch=nothing`: override the architecture (usually inferred from `grid`)
- `sinking_tracers=nothing`: sinking speed overrides, e.g. `(POC = 10/day, ...)`
- `open_bottom=true`: whether sinking tracers leave the domain

Returns
-------
`construct` returns an `Oceananigans.Biogeochemistry.AbstractContinuousFormBiogeochemistry` instance.
`construct_with_manifest` returns `(bgc, manifest)`, where `manifest` is a JSON-compatible construction record.

Example
-------
```julia
using Agate.Models: DARWIN

bgc = DARWIN.construct()
```

Trait-style allometric parameter overrides may be supplied during construction:

```julia
using Oceananigans.Units: day
using Agate.Library.Allometry: AllometricParam, PowerLaw

bgc = DARWIN.construct(;
    parameters=(;
        maximum_growth_rate=AllometricParam(
            PowerLaw(); prefactor=2 / day, exponent=-0.15
        ),
    ),
)
```
"""
function construct(; kwargs...)
    inputs = _construction_inputs(; kwargs...)
    return Construction.construct_factory(inputs.factory; inputs.kwargs...)
end

"""
    construct_with_manifest(; kw...) -> bgc, manifest

Construct a model instance and return it with a JSON-compatible construction manifest.
"""
function construct_with_manifest(; kwargs...)
    inputs = _construction_inputs(; kwargs...)
    bgc, resolved = Construction.construct_factory_with_context(
        inputs.factory; inputs.kwargs...
    )
    manifest = finalize_construction_manifest(
        default_model_manifest_context(
            :DARWIN, resolved; recipe=default_model_recipe(:DARWIN, resolved)
        )
    )
    return bgc, manifest
end
