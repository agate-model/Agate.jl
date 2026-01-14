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

# NOTE: This submodule lives under `Agate` (i.e. `Agate.Constructor`).
# Refer to sibling modules via relative paths (`..`) instead of the parent
# module name (`Agate`) because the parent name is not guaranteed to be bound
# inside submodules (especially for code that is `eval`'d into modules).
using ..Utils.Specifications: PFTSpecification, ModelSpecification, pft_get, pft_has

export construct
# Construction-time parameter registry bundles
export default_parameter_args, ParameterRegistryArgs


# Convenience update helpers
export update_plankton_args, update_dynamics

# Low-level patching is still available for advanced use
export patch

# Re-export key parameter containers as part of the constructor surface.
export PFTSpecification, ModelSpecification

include("patch.jl")
include("construct.jl")

end # module
