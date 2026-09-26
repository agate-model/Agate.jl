"""Integration layers between compiled Agate runtimes and external biogeochemistry frameworks."""
module Integrations

include("oceanbiome_npd.jl")

export NPDPlankton
export construct_npd_plankton, construct_npd_plankton_plus_recipe
export npd_configuration

end # module
