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
using ..Configuration: PFTSpecification

export construct_factory
export PFTSpecification

include("generator.jl")
include("construct.jl")

end # module
