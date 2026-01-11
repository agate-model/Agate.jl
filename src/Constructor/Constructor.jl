"""
Model-agnostic biogeochemistry constructor.

The public entry point is:

```julia
Agate.Constructor.construct(factory::AbstractBGCFactory; kwargs...) -> bgc_type
```

`bgc_type()` instantiates the biogeochemistry object with the constructed defaults.

This module also provides convenience helpers for updating the default argument
containers returned by `Agate.Models.*` factories.
"""
module Constructor

using Agate.Utils.Specifications: PFTSpecification, BiogeochemistrySpecification, ModelSpecification, pft_get, pft_has, cast_pft, cast_spec

export construct

# Convenience update helpers
export update_plankton_args, update_biogeochem_args, update_dynamics, update_group

# Low-level patching is still available for advanced use
export patch

# Re-export key parameter containers as part of the constructor surface.
export PFTSpecification, BiogeochemistrySpecification, ModelSpecification

include("patch.jl")
include("construct.jl")

end # module
