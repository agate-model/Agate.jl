"""
Module to construct an instance of a size-structured NiPiZD model.
"""

module Constructor

using Agate.Utils
using Agate.Models.Parameters
using Agate.Models.NiPiZD.Tracers

using UUIDs
using OceanBioME
using Oceananigans.Units

using OceanBioME: setup_velocity_fields
using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry

export construct, instantiate
export DEFAULT_PHYTO_ARGS,
    DEFAULT_PHYTO_GEIDER_ARGS, DEFAULT_ZOO_ARGS, DEFAULT_INTERACTION_ARGS, DEFAULT_BGC_ARGS

DEFAULT_PHYTO_ARGS = Dict(
    "allometry" => Dict(
        "maximum_growth_rate" => Dict("a" => 2 / day, "b" => -0.15),
        "nutrient_half_saturation" => Dict("a" => 0.17, "b" => 0.27),
    ),
    "linear_mortality" => 8e-7 / second,
    "alpha" => 0.1953 / day,
)

DEFAULT_PHYTO_GEIDER_ARGS = Dict(
    "allometry" => Dict(
        "maximum_growth_rate" => Dict("a" => 2 / day, "b" => -0.15),
        "nutrient_half_saturation" => Dict("a" => 0.17, "b" => 0.27),
    ),
    "linear_mortality" => 8e-7 / second,
    "photosynthetic_slope" => 0.46e-5,
    "chlorophyll_to_carbon_ratio" => 0.1,
)

DEFAULT_ZOO_ARGS = Dict(
    "allometry" => Dict("maximum_predation_rate" => Dict("a" => 30.84 / day, "b" => -0.16)),
    "linear_mortality" => 8e-7 / second,
    "holling_half_saturation" => 5.0,
    "quadratic_mortality" => 1e-6 / second,
)

DEFAULT_INTERACTION_ARGS = Dict(
    "P" => Dict(
        "can_eat" => 0, # bool
        "can_be_eaten" => 1, # bool
        "optimum_predator_prey_ratio" => 0,
        "protection" => 0,
        "specificity" => 0,
        "assimilation_efficiency" => 0,
    ),
    "Z" => Dict(
        "can_eat" => 1, # bool
        "can_be_eaten" => 0, # bool
        "optimum_predator_prey_ratio" => 10,
        "protection" => 1,
        "specificity" => 0.3,
        "assimilation_efficiency" => 0.32,
    ),
)

DEFAULT_BGC_ARGS = Dict(
    "detritus_remineralization" => 0.1213 / day, "mortality_export_fraction" => 0.5
)

"""
    construct(;
        n_phyto=2,
        n_zoo=2,
        phyto_diameters=Dict(
            "min_diameter" => 2, "max_diameter" => 10, "splitting" => "log_splitting"
        ),
        zoo_diameters=Dict(
            "min_diameter" => 20, "max_diameter" => 100, "splitting" => "linear_splitting"
        ),
        nutrient_dynamics=nutrients_default,
        detritus_dynamics=detritus_default,
        phyto_dynamics=phytoplankton_default,
        zoo_dynamics=zooplankton_default,
        phyto_args=DEFAULT_PHYTO_ARGS,
        zoo_args=DEFAULT_ZOO_ARGS,
        interaction_args=DEFAULT_INTERACTION_ARGS,
        bgc_args=DEFAULT_BGC_ARGS,
        palatability_matrix=nothing,
        assimilation_efficiency_matrix=nothing,
        sinking_tracers=nothing,
        grid=BoxModelGrid(),
        open_bottom=true,
    )

Construct a size-structured NiPiZD model abstract type.

This constructor builds a size-structured plankton model with two plankton functional types:
phytoplankton (P) and zooplankton (Z), each of which can be specified to have any number of
size classes (`n_phyto` and `n_zoo`). In addition to plankton, the constructor implements
idealized detritus (D) and nutrient (N) cycling by default, although more complex N and D
cycling can also be defined using the `nutrient_dynamics` and `detritus_dynamics` arguments.

During model construction, the size of each plankton determines photosynthetic growth rates,
nutrient half saturation constants, predation rates, and predator-prey assimilation and
palatability values. Alternatively, if manually defined predator-prey assimilation and
palatability values are desired, these can be specified using the `palatability_matrix` and
`assimilation_efficiency_matrix` arguments.

Note that if non-default `*_dynamics` expressions are passed, the relevant `*_args` also
need to be specified.

The type specification includes a photosynthetic active radiation (PAR) auxiliary field.

# Keywords
- `n_phyto`: number of phytoplankton in the model
- `n_zoo`: number of zooplankton in the model
- `phyto_diameters`: dictionary from which `phyto` diameters can be computed or a list of
    values to use (as many as the model expects)
- `zoo_diameters`: dictionary from which `zoo` diameters can be computed or a list of
    values to use (as many as the model expects)
- `nutrient_dynamics`: expression describing how nutrients change over time, see
    `Agate.Models.Tracers`
- `detritus_dynamics`: expression describing how detritus evolves over time, see
    `Agate.Models.Tracers`
- `phyto_dynamics`: expression describing how phytoplankton grow, see `Agate.Models.Tracers`
- `zoo_dynamics`: expression describing how zooplankton grow, see `Agate.Models.Tracers`
- `phyto_args`: Dictionary of phytoplankton parameters, for default values see
    `NiPiZD.DEFAULT_PHYTO_ARGS`
- `zoo_args`: Dictionary of zooplankton parameters, for default values see
    `NiPiZD.DEFAULT_ZOO_ARGS`
- `interaction_args`: Dictionary of arguments from which a palatability and assimilation
   efficiency matrix between all plankton can be computed, for default values see
    `NiPiZD.DEFAULT_INTERACTION_ARGS`
- `bgc_args`: Dictionary of biogeochemistry parameters related to nutrient and detritus, for
    default values see `NiPiZD.DEFAULT_BGC_ARGS`
- `palatability_matrix`: optional palatability matrix passed as an Array, if provided
    then `interaction_args` are not used to compute this
- `assimilation_efficiency_matrix`: optional assimilation efficiency matrix passed as an
    Array, if provided then `interaction_args` are not used to compute this
- `sinking_tracers`: optional NamedTuple of sinking speeds (passed as positive values) of
   the form (<tracer name expressed as symbol> = <speed>, ...)
- `grid`: optional Oceananigans grid object defining the geometry to build the model on, must
   be passed if `sinking_tracers` is defined, defaults to BoxModelGrid
- `open_bottom`: indicates whether the sinking velocity should be smoothly brought to zero
   at the bottom to prevent the tracers leaving the domain, defaults to `true`, which means
   the bottom is open and the tracers leave (i.e., no slowing of velocity to 0 is applied)

# Example
```julia
using Agate.Constructors: NiPiZD

n2p2zd = NiPiZD.construct()
n2p2zd_model_obj = n2p2zd()
```
"""
function construct(;
    n_phyto=2,
    n_zoo=2,
    phyto_diameters=Dict(
        "min_diameter" => 2, "max_diameter" => 10, "splitting" => "log_splitting"
    ),
    zoo_diameters=Dict(
        "min_diameter" => 20, "max_diameter" => 100, "splitting" => "linear_splitting"
    ),
    nutrient_dynamics=nutrients_default,
    detritus_dynamics=detritus_default,
    phyto_dynamics=phytoplankton_default,
    zoo_dynamics=zooplankton_default,
    phyto_args=DEFAULT_PHYTO_ARGS,
    zoo_args=DEFAULT_ZOO_ARGS,
    interaction_args=DEFAULT_INTERACTION_ARGS,
    bgc_args=DEFAULT_BGC_ARGS,
    palatability_matrix=nothing,
    assimilation_efficiency_matrix=nothing,
    sinking_tracers=nothing,
    grid=BoxModelGrid(),
    open_bottom=true,
)
    parameters, plankton_names = create_size_structured_params(;
        n_plankton=Dict("P" => n_phyto, "Z" => n_zoo),
        diameters=Dict("P" => phyto_diameters, "Z" => zoo_diameters),
        plankton_args=Dict("P" => phyto_args, "Z" => zoo_args),
        interaction_args=interaction_args,
        bgc_args=bgc_args,
        palatability_matrix=palatability_matrix,
        assimilation_efficiency_matrix=assimilation_efficiency_matrix,
    )

    # NOTE: Zs precede Ps
    plankton_array = [Symbol(name) for name in plankton_names]

    # create tracer functions
    tracers = Dict(
        "N" => nutrient_dynamics(plankton_array), "D" => detritus_dynamics(plankton_array)
    )

    for i in 1:n_zoo
        name = "Z$i"
        index = findfirst(x -> x == name, plankton_names)
        tracers[name] = zoo_dynamics(plankton_array, name, index)
    end

    for i in 1:n_phyto
        name = "P$i"
        index = findfirst(x -> x == name, plankton_names)
        tracers[name] = phyto_dynamics(plankton_array, name, index)
    end

    if !isnothing(sinking_tracers)
        sinking_velocities = setup_velocity_fields(sinking_tracers, grid, open_bottom)
    else
        sinking_velocities = nothing
    end

    # return Oceananigans.Biogeochemistry object
    # NOTE: this adds "PAR" as an auxiliary field by default
    return define_tracer_functions(
        parameters, tracers; sinking_velocities=sinking_velocities
    )
end

"""
    instantiate(
        bgc_type;
        phyto_diameters=Dict(
            "min_diameter" => 2, "max_diameter" => 10, "splitting" => "log_splitting"
        ),
        zoo_diameters=Dict(
            "min_diameter" => 20, "max_diameter" => 100, "splitting" => "linear_splitting"
        ),
        phyto_args=DEFAULT_PHYTO_ARGS,
        zoo_args=DEFAULT_ZOO_ARGS,
        interaction_args=DEFAULT_INTERACTION_ARGS,
        bgc_args=DEFAULT_BGC_ARGS,
        palatability_matrix=nothing,
        assimilation_efficiency_matrix=nothing,
        sinking_tracers=nothing,
        grid=BoxModelGrid(),
        open_bottom=true,
    )

A function to instantiate an object of `bgc_type` returned by `NiPiZD.construct()`.

The type specifies the number of phytoplankton and zooplankton in the model and includes
default parameter values. The instantiate method is used to override the default values
of any of the model parameters or plankton diameters.

# Arguments
- `bgc_type`: subtype of Oceananigans.Biogeochemistry returned by `NiPiZD.construct()`
   with a specified number of phytoplankton and zooplankton

# Keywords
- `phyto_diameters`: dictionary from which `phyto` diameters can be computed or a list of
    values to use (as many as the model expects)
- `zoo_diameters`: dictionary from which `zoo` diameters can be computed or a list of
    values to use (as many as the model expects)
- `nutrient_dynamics`: expression describing how nutrients change over time, see
    `Agate.Models.Tracers`
- `detritus_dynamics`: expression describing how detritus evolves over time, see
    `Agate.Models.Tracers`
- `phyto_dynamics`: expression describing how phytoplankton grow, see `Agate.Models.Tracers`
- `zoo_dynamics`: expression describing how zooplankton grow, see `Agate.Models.Tracers`
- `phyto_args`: Dictionary of phytoplankton parameters, for default values see
    `NiPiZD.DEFAULT_PHYTO_ARGS`
- `zoo_args`: Dictionary of zooplankton parameters, for default values see
    `NiPiZD.DEFAULT_ZOO_ARGS`
- `interaction_args`: Dictionary of arguments from which a palatability and assimilation
   efficiency matrix between all plankton can be computed, for default values see
    `NiPiZD.DEFAULT_INTERACTION_ARGS`
- `bgc_args`: Dictionary of constant parameters used in growth functions (i.e., not size
    dependant plankton parameters as well as biogeochemistry parameters related to nutrient
    and detritus, for default values see `NiPiZD.DEFAULT_CONSTANT_ARGS`
- `palatability_matrix`: optional palatability matrix passed as an Array, if provided
    then `interaction_args` are not used to compute this
- `assimilation_efficiency_matrix`: optional assimilation efficiency matrix passed as an
    Array, if provided then `interaction_args` are not used to compute this
- `sinking_tracers`: optional NamedTuple of sinking speeds (passed as positive values) of
   the form (<tracer name expressed as symbol> = <speed>, ...)
- `grid`: optional Oceananigans grid object defining the geometry to build the model on, must
   be passed if `sinking_tracers` is defined, defaults to BoxModelGrid
- `open_bottom`: indicates whether the sinking velocity should be smoothly brought to zero
   at the bottom to prevent the tracers leaving the domain, defaults to `true`, which means
   the bottom is open and the tracers leave (i.e., no slowing of velocity to 0 is applied)

# Example
```julia
using Agate.Constructors: NiPiZD

n2p2zd = NiPiZD.construct()

# change some parameter values
phyto_args = deepcopy(NiPiZD.DEFAULT_PHYTO_ARGS)
phyto_args["allometry"]["maximum_growth_rate"]["a"] = 2
n2p2zd_model_obj = NiPiZD.instantiate(n2p2zd; phyto_args=phyto_args)
```
"""
function instantiate(
    bgc_type;
    phyto_diameters=Dict(
        "min_diameter" => 2, "max_diameter" => 10, "splitting" => "log_splitting"
    ),
    zoo_diameters=Dict(
        "min_diameter" => 20, "max_diameter" => 100, "splitting" => "linear_splitting"
    ),
    phyto_args=DEFAULT_PHYTO_ARGS,
    zoo_args=DEFAULT_ZOO_ARGS,
    interaction_args=DEFAULT_INTERACTION_ARGS,
    bgc_args=DEFAULT_BGC_ARGS,
    palatability_matrix=nothing,
    assimilation_efficiency_matrix=nothing,
    sinking_tracers=nothing,
    grid=BoxModelGrid(),
    open_bottom=true,
)
    defaults = bgc_type()
    n_phyto = Int(defaults.n_P)
    n_zoo = Int(defaults.n_Z)

    parameters, _ = create_size_structured_params(;
        n_plankton=Dict("P" => n_phyto, "Z" => n_zoo),
        diameters=Dict("P" => phyto_diameters, "Z" => zoo_diameters),
        plankton_args=Dict("P" => phyto_args, "Z" => zoo_args),
        interaction_args=interaction_args,
        bgc_args=bgc_args,
        palatability_matrix=palatability_matrix,
        assimilation_efficiency_matrix=assimilation_efficiency_matrix,
    )

    # params are a NamedTuple -> have to convert to Dict
    parameters_dict = Dict(pairs(parameters))

    if !isnothing(sinking_tracers)
        if !(:sinking_velocities in fieldnames(bgc_type))
            throw(ArgumentError("BGC object does not have sinking_velocities"))
        end
        if isnothing(grid)
            throw(ArgumentError("grid must be defined to setup tracer sinking"))
        end
        parameters_dict[:sinking_velocities] = setup_velocity_fields(
            sinking_tracers, grid, open_bottom
        )
    end

    return bgc_type(; parameters_dict...)
end

end # module
