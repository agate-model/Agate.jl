"""
Model-agnostic biogeochemistry constructor.

The public entry point is:

```julia
Agate.Constructor.construct(factory::AbstractBGCFactory; kwargs...) -> bgc
```

`construct` returns the biogeochemistry instance directly.

This module also provides convenience helpers for updating the default argument
containers returned by `Agate.Models.*` factories.
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
## Convenience update helpers
export update_community, extend_community
export update_dynamics, extend_dynamics

# Low-level patching is still available for advanced use
export patch

# Re-export key parameter containers as part of the constructor surface.
export PFTSpecification, ModelSpecification

include("model_spec.jl")
include("patch.jl")
include("generator.jl")
include("construct.jl")

end # module
