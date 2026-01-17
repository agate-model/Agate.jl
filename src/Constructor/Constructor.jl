"""
Model-agnostic biogeochemistry constructor.

The public entry points are the **model constructors** exposed in `Agate.Models`
(e.g. `Agate.Models.NiPiZD.NiPiZD(...)` and `Agate.Models.DARWIN.DARWIN(...)`).

This module contains the model-agnostic implementation:

```julia
Agate.Constructor.construct(factory::AbstractBGCFactory; kwargs...) -> bgc
```

`construct` returns the biogeochemistry instance directly.

Design goals
------------
- Maintainability over cleverness.
- Explicit behavior (no hidden magic).
- GPU compatibility via `Adapt.jl` + Oceananigans architectures.
- Preserve composability with Oceananigans/OceanBioME conventions.

Anything that *patches* or *extends* configuration containers is intentionally
left to normal Julia mechanisms (`merge`, `NamedTuple` construction, etc.) to keep
the public surface area small.
"""
module Constructor

# Oceananigans is a core dependency for Agate's host/GPU architecture model.
# Import it here so constructor code can reliably query `Oceananigans.Architectures`.
import Oceananigans

# NOTE: This submodule lives under `Agate` (i.e. `Agate.Constructor`).
# Refer to sibling modules via relative paths (`..`) instead of the parent
# module name (`Agate`) because the parent name is not guaranteed to be bound
# inside submodules (especially for code that is `eval`'d into modules).
using ..Utils.Specifications: PFTSpecification, ModelSpecification

export construct

# Re-export key parameter containers as part of the constructor surface.
export PFTSpecification, ModelSpecification

include("model_spec.jl")
include("generator.jl")
include("construct.jl")

end # module
