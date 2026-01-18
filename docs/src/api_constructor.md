# Constructor API

Agate's primary user-facing constructors are model-specific and live under the
model modules:

- `NiPiZD.construct`
- `DARWIN.construct`

Each constructor returns a concrete biogeochemistry **instance** (already
`Adapt.jl` compatible) that composes with OceanBioME/Oceananigans.

## Core workflow

```julia
using Agate

bgc = NiPiZD.construct()  # CPU instance by default
```

Precision is determined by the **grid element type** when you pass `grid=...`.
When no grid is provided (or when using the default `BoxModelGrid()`), the instance
uses `Float64`.

To construct directly on a GPU architecture, either:

- request a GPU instance via `arch=GPU(...)`, or
- adapt the returned instance explicitly (for example with CUDA).

```julia
using Adapt
using CUDA
using Oceananigans.Architectures: GPU

bgc = NiPiZD.construct()
bgc_gpu = Adapt.adapt(CuArray, bgc)

# Or request a GPU instance directly.
bgc_gpu2 = NiPiZD.construct(; arch=GPU())
```

When interfacing with Oceananigans/OceanBioME, it is usually better to let the
host grid/model choose precision and architecture. Pass the grid you will use;
Agate infers precision from `eltype(grid)` and defaults `arch = architecture(grid)`.

```julia
using Oceananigans
using Oceananigans.Architectures: GPU

grid = RectilinearGrid(size=(16, 16, 16), extent=(1, 1, 1),
                       topology=(Periodic, Periodic, Bounded),
                       architecture=GPU(), float_type=Float32)

bgc_gpu = NiPiZD.construct(; grid)  # precision + architecture inferred from `grid`
```

## Parameter overrides

Model constructors accept `parameters=(; ...)` as a **NamedTuple of overrides**.
Unknown parameter keys throw immediately (typo protection).

```julia
using Oceananigans.Units: day

bgc = NiPiZD.construct(; parameters=(detritus_remineralization = 0.18 / day,))
```

## NiPiZD interaction overrides

NiPiZD allows overriding the two interaction matrices directly.

You can pass **concrete matrices**:

```julia
using Agate

pal = [0.0 0.0 1.0 0.25;
       0.0 0.0 0.10 1.0]        # n_zoo x n_phyto

assim = fill(0.32, size(pal))

bgc = NiPiZD.construct(; palatability_matrix=pal, assimilation_matrix=assim)
```

Or, for a beginner-friendly workflow, pass **provider functions** that take only
diameters and group names (no `ctx ->` required). One ergonomic pattern is to
capture `n_phyto`/`n_zoo` in a closure:

```julia
n_phyto = 4
n_zoo = 2

my_palatability = function (diameters, group_symbols)
    pal = zeros(Float32, n_zoo, n_phyto)
    # ... fill palatability here ...
    return pal
end

my_assimilation = function (diameters, group_symbols)
    return fill(0.32f0, n_zoo, n_phyto)
end

bgc = NiPiZD.construct(; n_phyto, n_zoo,
                         palatability_matrix=my_palatability,
                         assimilation_matrix=my_assimilation)
```

## Introspection helpers

For interactive work (REPL/notebooks), Agate provides small helpers to discover
the model surface:

```julia
names = tracer_names(bgc)                 # Vector{Symbol}
aux   = auxiliary_field_names(bgc)        # Vector{Symbol}
pars  = parameter_names(bgc)              # Vector{Symbol}

describe(bgc)                            # prints a short summary
summary = model_summary(bgc)              # returns a NamedTuple
```

## API reference

```@docs
Agate.NiPiZD.construct
Agate.DARWIN.construct
Agate.tracer_names
Agate.auxiliary_field_names
Agate.model_summary
Agate.describe
Agate.parameter_names
Agate.required_parameters
```
