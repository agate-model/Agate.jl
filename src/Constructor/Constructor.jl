"""Model-agnostic biogeochemistry constructor.

The public entry point is:

```julia
Agate.Constructor.construct(factory::AbstractBGCFactory; kwargs...) -> bgc_type
```

`bgc_type()` instantiates the biogeochemistry object with the constructed defaults.

Utilities for ergonomically overriding factory configuration.

Agate uses immutable configuration objects (primarily `NamedTuple`s and lightweight
containers around `NamedTuple`s). This file provides a small, non-mutating
"patch" layer that creates overridden copies and (for `NamedTuple`s) checks that
keys exist to prevent silent typos.
"""
module Constructor

using Agate.Utils.Specifications: PFTSpecification, BiogeochemistrySpecification, ModelSpecification, pft_get, pft_has, cast_pft, cast_spec

export construct
export patch, update_group

# Re-export key parameter containers as part of the constructor surface.
export PFTSpecification, BiogeochemistrySpecification, ModelSpecification

include("patch.jl")
include("construct.jl")

end # module
