"""
Model-agnostic biogeochemistry constructor.

The public entry points are the **model constructors** exposed in `Agate.Models`
(e.g. `Agate.Models.NiPiZD.NiPiZD(...)` and `Agate.Models.DARWIN.DARWIN(...)`).

This module contains the model-agnostic implementation:

```julia
Agate.Construction.construct_factory(factory::AbstractBGCFactory; kwargs...) -> bgc
```

`construct_factory` returns the biogeochemistry instance directly.

"""
module Construction

import Oceananigans

# NOTE: This submodule lives under `Agate` (i.e. `Agate.Construction`).
# Refer to sibling modules via relative paths (`..`) instead of the parent
# module name (`Agate`) because the parent name is not guaranteed to be bound
# inside submodules (especially for code that is `eval`'d into modules).
using ..Configuration: PFTSpecification

export construct_factory

# Re-export key parameter containers as part of the constructor surface.
export PFTSpecification
include("generator.jl")
include("construct.jl")

end # module
