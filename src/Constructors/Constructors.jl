"""
Module for high level model constructors.
"""

module Constructors

include("NPZD_size_structured.jl")
include("thunder_egg_1.jl")

using .NPZD_size_structured
using .thunder_egg_1

export construct_size_structured_NPZD
export construct_thunder_egg_1

end # module
