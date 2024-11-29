module Constructors

export construct_NPZD_instance

function construct_NPZD_instance(
    n_phyto=1,
    n_zoo=1,
    # nutrient_tracers ??
    detritus=typical_detritus,
    phyto_growth=phytoplankton_growth,
    zoo_growth=zooplankton_growth,
    phyto_args=Dict(),
    zoo_args=Dict(),
    palatability_matrix=nothing,
    assimilation_efficiency_matrix=nothing,
) end

end # module
