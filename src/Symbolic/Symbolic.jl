"""Construction-time symbolic equation builders.

These helpers produce `Agate.Equations.Equation` / `Agate.Equations.AExpr` objects.
They are intended for model authoring; generated kernels operate on numeric state.
"""

module Symbolic

include("Mortality.jl")
include("Predation.jl")
include("Photosynthesis.jl")
include("Remineralization.jl")

using .Mortality
using .Predation
using .Photosynthesis
using .Remineralization

end # module
