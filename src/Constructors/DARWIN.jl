"""
Module to construct an instance of an Agate.jl-DARWIN model.
"""

module DARWIN

using Agate.Models.Biogeochemistry
using Agate.Models.Parameters
using Agate.Models.Tracers

using NamedArrays
using UUIDs
using Oceananigans.Units

using OceanBioME: setup_velocity_fields
using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry

DEFAULT_PHYTO_ARGS = Dict(
    "allometry" => Dict(
        "maximum_growth_rate" => Dict("a" => 2 / day, "b" => -0.15),
        "half_saturation_DIN" => Dict("a" => 0.17, "b" => 0.27),
        "half_saturation_PO4" => Dict("a" => 0.17, "b" => 0.27),
        # need this to vectorize the tracer functions
        "maximum_predation_rate" => Dict("a" => 0, "b" => 0),
    ),
    "linear_mortality_p" => 8e-7 / second,
    "photosynthetic_slope" => 0.46e-5,
    "chlorophyll_to_carbon_ratio" => 0.1,
)

DEFAULT_ZOO_ARGS = Dict(
    "allometry" => Dict("maximum_predation_rate" => Dict("a" => 30.84 / day, "b" => -0.16)),
    "linear_mortality_z" => 8e-7 / second,
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
    "POC_remineralization" => 0.1213 / day,
    "DOC_remineralization" => 0.1213 / day,
    "PON_remineralization" => 0.1213 / day,
    "DON_remineralization" => 0.1213 / day,
    "POP_remineralization" => 0.1213 / day,
    "DOP_remineralization" => 0.1213 / day,
    "DOM_POM_fractionation" => 0.45,
    "nitrogen_to_carbon" => 0.15,
    "phosphorus_to_carbon" => 0.009,
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
        DIC_dynamics=DIC_geider_light,
        PO4_dynamics=PO4_geider_light_fixed_ratios,
        DIN_dynamics=DIN_geider_light_fixed_ratios,
        POC_dynamics=POC_typical,
        DOC_dynamics=DOC_typical,
        PON_dynamics=PON_typical,
        DON_dynamics=DON_typical,
        POP_dynamics=POP_typical,
        DOP_dynamics=DOP_typical,
        phyto_dynamics=phytoplankton_growth_two_nutrients_geider_light,
        zoo_dynamics=zooplankton_growth_simplified,
        phyto_args=DEFAULT_PHYTO_ARGS,
        zoo_args=DEFAULT_ZOO_ARGS,
        interaction_args=DEFAULT_INTERACTION_ARGS,
        bgc_args=DEFAULT_BGC_ARGS,
        palatability_matrix=nothing,
        assimilation_efficiency_matrix=nothing,
    ) -> DataType

Construct an Agate.jl-DARWIN model abstract type.

!!! info

    This model is in active development and has not been validated against `MITgcm-DARWIN`.

!!! formulation

    TRACERS:

    ∂t cⱼ = ``Uⱼ``DIC - ``Mⱼ`` + ``Gⱼ`` - ``gⱼ``

    ∂t DIC = ∑(``Uⱼ`` DIC) + ``R``DOC + ``R``POC

    ∂t DIN = ∑(``Uⱼ``DIC * ``Qⱼ``N)  + ``R``DON + ``R``PON

    ∂t PO4 = ∑(``Uⱼ``DIC * ``Qⱼ``P)  + ``R``DOP + ``R``POP

    ∂t DOC = ∑(``Mⱼ``DOC) + ``g``DOC - ``R``DOC

    ∂t DON = ∑(``Mⱼ``DOC * ``Qⱼ``N) + ``g``DON - ``R``DON

    ∂t DOP = ∑(``Mⱼ``DOC * ``Qⱼ``P) + ``g``DOP - ``R``DOP

    ∂t POC = ∑(``Mⱼ``POC) + ``g``POC - ``R``POC

    ∂t PON = ∑(``Mⱼ``POC * ``Qⱼ``N) + ``g``PON - ``R``PON

    ∂t POP = ∑(``Mⱼ``POC * ``Qⱼ``P) + ``g``POP - ``R``POP

    where:
    - ``U`` = uptake
    - ``R`` = remineralization
    - ``M`` = mortality
    - ``g, G`` = grazing losses and gains
    - ``Q`` = plankton elemental ratios

    TRAITS:

    μmax, KR, gmax = a*Volume^b

    palat = η/(1+(``ratio``-``opt``)^2)^σ

    where:
    - μmax = maximum photosynthetic growth
    - KR = nutrient half saturation
    - gmax = maximum predation rate
    - palat = palatability
    - ``ratio`` = predator to prey size ratio (diameter)
    - ``opt`` = predator to prey size optimum (diameter)
    - η = prey protection
    - σ = predator specificity


This constructor builds an Agate.jl-DARWIN model with two plankton functional types:
phytoplankton (P) and zooplankton (Z), each of which can be specified to have any number of
size classes (`n_phyto` and `n_zoo`). In addition to plankton, the constructor implements
idealized dissolved inorganic carbon (DIC), particulate organic matter (POC, POP, PON), dissolved organic
matter (DOC, DOP, DON) and two nutrients (DIN and PO4) cycling by default, although more complex elemental cycling can also be defined using the `nutrient_dynamics` and `detritus_dynamics`
arguments.

During model construction, the size of each plankton determines photosynthetic growth rates,
nutrient half saturation constants, predation rates, and optionally predator-prey assimilation
and palatability values. Alternatively, if manually defined predator-prey assimilation and
palatability values are desired, these can be defined using the `palatability_matrix` and
`assimilation_efficiency_matrix` arguments.

Note that if non-default `*_dynamics` expressions are passed, the relevant `*_args` also
need to be specified.

The type specification includes a photosynthetic active radiation (PAR) auxiliary field.

# Arguments
- `n_phyto`: number of phytoplankton to include in the model
- `n_zoo`: number of zooplankton to include in the model
- `phyto_diameters`: dictionary from which `n_phyto` diameters can be computed or a list of
    values to use
- `zoo_diameters`: dictionary from which `zoo` diameters can be computed or a list of
    values to use
- `DIC_dynamics`: expression describing how DIC changes over time, see `Agate.Models.Tracers`
- `PO4_dynamics`: expression describing how PO4 changes over time, see `Agate.Models.Tracers`
- `DIN_dynamics`: expression describing how DIN changes over time, see `Agate.Models.Tracers`
- `POC_dynamics`: expression describing how POC changes over time, see `Agate.Models.Tracers`
- `DOC_dynamics`: expression describing how DOC changes over time, see `Agate.Models.Tracers`
- `PON_dynamics`: expression describing how PON changes over time, see `Agate.Models.Tracers`
- `DON_dynamics`: expression describing how DON changes over time, see `Agate.Models.Tracers`
- `POP_dynamics`: expression describing how POP changes over time, see `Agate.Models.Tracers`
- `DOP_dynamics`: expression describing how DOP changes over time, see `Agate.Models.Tracers`
- `phyto_dynamics`: expression describing how phytoplankton grow, see `Agate.Models.Tracers`
- `zoo_dynamics`: expression describing how zooplankton grow, see `Agate.Models.Tracers`
- `phyto_args`: Dictionary of phytoplankton parameters, for default values see
    `Agate.Models.Constructors.DEFAULT_PHYTO_ARGS`
- `zoo_args`: Dictionary of zooplankton parameters, for default values see
    `Agate.Models.Constructors.DEFAULT_ZOO_ARGS`
- `interaction_args`: Dictionary of arguments from which a palatability and assimilation
   efficiency matrix between all plankton can be computed, for default values see
    `Agate.Models.Constructors.DEFAULT_INTERACTION_ARGS`
- `bgc_args`: biogeochemistry parameters related to nutrient and detritus, for default
    values see `Agate.Models.Constructors.DEFAULT_BGC_ARGS`
- `palatability_matrix`: optional palatability matrix passed as a NamedArray, if provided
   then `interaction_args` are not used to compute this
- `assimilation_efficiency_matrix`: optional assimilation efficiency matrix passed as a
   NamedArray, if provided then `interaction_args` are not used to compute this


# Example
```julia
using Agate.Constructors: DARWIN

darwin_2p_2z = DARWIN.construct()
darwin_2p_2z_model_obj = darwin_2p_2z()
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
    DIC_dynamics=DIC_geider_light,
    PO4_dynamics=PO4_geider_light_fixed_ratios,
    DIN_dynamics=DIN_geider_light_fixed_ratios,
    POC_dynamics=POC_typical,
    DOC_dynamics=DOC_typical,
    PON_dynamics=PON_typical,
    DON_dynamics=DON_typical,
    POP_dynamics=POP_typical,
    DOP_dynamics=DOP_typical,
    phyto_dynamics=phytoplankton_growth_two_nutrients_geider_light,
    zoo_dynamics=zooplankton_growth_simplified,
    phyto_args=DEFAULT_PHYTO_ARGS,
    zoo_args=DEFAULT_ZOO_ARGS,
    interaction_args=DEFAULT_INTERACTION_ARGS,
    bgc_args=DEFAULT_BGC_ARGS,
    palatability_matrix=nothing,
    assimilation_efficiency_matrix=nothing,
)
    parameters = create_params_dict(;
        n_phyto=n_phyto,
        n_zoo=n_zoo,
        phyto_diameters=phyto_diameters,
        zoo_diameters=zoo_diameters,
        phyto_args=phyto_args,
        zoo_args=zoo_args,
        interaction_args=interaction_args,
        bgc_args=bgc_args,
        palatability_matrix=palatability_matrix,
        assimilation_efficiency_matrix=assimilation_efficiency_matrix,
    )

    # create tracer functions
    plankton_array = vcat(
        [Symbol("Z$i") for i in 1:n_zoo], [Symbol("P$i") for i in 1:n_phyto]
    )
    tracers = Dict(
        "DIC" => DIC_dynamics(plankton_array),
        "DIN" => DIN_dynamics(plankton_array),
        "PO4" => PO4_dynamics(plankton_array),
        "POC" => POC_dynamics(plankton_array),
        "DOC" => DOC_dynamics(plankton_array),
        "PON" => PON_dynamics(plankton_array),
        "DON" => DON_dynamics(plankton_array),
        "POP" => POP_dynamics(plankton_array),
        "DOP" => DOP_dynamics(plankton_array),
    )
    for i in 1:n_phyto
        name = "P$i"
        tracers[name] = phyto_dynamics(plankton_array, name)
    end
    for i in 1:n_zoo
        name = "Z$i"
        tracers[name] = zoo_dynamics(plankton_array, name)
    end

    # return Oceananigans.Biogeochemistry object
    # note this adds "PAR" as an auxiliary field by default
    return define_tracer_functions(parameters, tracers)
end

"""
    instantiate(
        bgc_type;
        n_phyto=2,
        n_zoo=2,
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
    ) -> bgc_type

A function to instantiate an object of `bgc_type` returned by `DARWIN.construct()`.

The type specifies the number of phytoplankton and zooplankton in the model and includes
default parameter values. The instantiate method is used to override the default values
of any of the model parameters or plankton diameters.

!!! tip

    Changing the parameter values of an existing DARWIN model type using `instantiate()` is useful in
    dynamic programming contexts such as `for` loops.

# Arguments
- `bgc_type`: subtype of Oceananigans.Biogeochemistry returned by `DARWIN.construct()`
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
    `Agate.Models.Constructors.DEFAULT_PHYTO_ARGS`
- `zoo_args`: Dictionary of zooplankton parameters, for default values see
    `Agate.Models.Constructors.DEFAULT_ZOO_ARGS`
- `interaction_args`: Dictionary of arguments from which a palatability and assimilation
   efficiency matrix between all plankton can be computed, for default values see
    `Agate.Models.Constructors.DEFAULT_INTERACTION_ARGS`
- `bgc_args`: Dictionary of constant parameters used in growth functions (i.e., not size
    dependant plankton parameters as well as biogeochemistry parameters related to nutrient
    and detritus, for default values see `Agate.Models.Constructors.DEFAULT_CONSTANT_ARGS`
- `palatability_matrix`: optional palatability matrix passed as a NamedArray, if provided
    then `interaction_args` are not used to compute this
- `assimilation_efficiency_matrix`: optional assimilation efficiency matrix passed as a
    NamedArray, if provided then `interaction_args` are not used to compute this

# Example
```julia
using Agate.Constructors: DARWIN

darwin_2p_2z = DARWIN.construct()

# change some parameter values
phyto_args = DARWIN.DEFAULT_PHYTO_ARGS
phyto_args["allometry"]["maximum_growth_rate"]["a"] = 2
darwin_2p_2z_model_obj = DARWIN.instantiate(darwin_2p_2z; phyto_args=phyto_args)
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
)
    defaults = bgc_type()
    n_phyto = Int(defaults.n_phyto)
    n_zoo = Int(defaults.n_zoo)

    # returns NamedTuple -> have to convert to Dict
    parameters = create_params_dict(;
        n_phyto=n_phyto,
        n_zoo=n_zoo,
        phyto_diameters=phyto_diameters,
        zoo_diameters=zoo_diameters,
        phyto_args=phyto_args,
        zoo_args=zoo_args,
        interaction_args=interaction_args,
        bgc_args=bgc_args,
        palatability_matrix=palatability_matrix,
        assimilation_efficiency_matrix=assimilation_efficiency_matrix,
    )

    return bgc_type(; Dict(pairs(parameters))...)
end

end # module
