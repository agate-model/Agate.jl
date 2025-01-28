"""
Module for high level model constructors.
"""

module Constructors

include("NPZD_size_structured.jl")
include("NPZD_two_nutrients.jl")

using .NPZD_size_structured
using .NPZD_two_nutrients

export construct_size_structured_NPZD
export construct_size_structured_P_Z_POC_DOC_DIN_PO4

end # module
