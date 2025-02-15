"""
Module to construct an instance of an Thunder Egg 1 model.
"""

module thunder_egg_1

using Agate.Models.Biogeochemistry
using Agate.Models.Parameters
using Agate.Models.Tracers

using Oceananigans.Units

export construct_thunder_egg_1

DEFAULT_PHYTO_ARGS = Dict(
    "allometry" => Dict(
        "maximum_growth_rate" => Dict("a" => 2 / day, "b" => -0.15),
        "half_saturation_DIN" => Dict("a" => 0.17, "b" => 0.27),
        "half_saturation_PO4" => Dict("a" => 0.17, "b" => 0.27),
    ),
    "linear_mortality" => 8e-7 / second,
    "photosynthetic_slope" => 0.46e-5,
    "chlorophyll_to_carbon_ratio" => 0.1,
    "nitrogen_to_carbon" => 0.15,
    "phosphorus_to_carbon" => 0.009,
)

DEFAULT_ZOO_ARGS = Dict(
    "allometry" => Dict("maximum_predation_rate" => Dict("a" => 30.84 / day, "b" => -0.16)),
    "linear_mortality" => 8e-7 / second,
    "holling_half_saturation" => 5.0,
    "quadratic_mortality" => 1e-6 / second,
    "nitrogen_to_carbon" => 0.15,
    "phosphorus_to_carbon" => 0.009,
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
        )

Construct an `Agate.jl-DARWIN` model object which is based on the `MITgcm-DARWIN` model.

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

    μmax, KN, gmax = a*Volume^b
    
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


This constructor builds a size-structured plankton model with two plankton functional types:
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
"""
function construct_thunder_egg_1(;
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
    phyto_args["n"] = n_phyto
    phyto_args["diameters"] = phyto_diameters
    zoo_args["n"] = n_zoo
    zoo_args["diameters"] = zoo_diameters

    # compute emergent parameters
    defined_parameters = Dict("P" => phyto_args, "Z" => zoo_args)

    if isnothing(palatability_matrix)
        defined_parameters["P"]["palatability"] = Dict(
            k => interaction_args["P"][k] for
            k in ["can_eat", "optimum_predator_prey_ratio", "protection", "specificity"]
        )
        defined_parameters["Z"]["palatability"] = Dict(
            k => interaction_args["Z"][k] for
            k in ["can_eat", "optimum_predator_prey_ratio", "protection", "specificity"]
        )
    end

    if isnothing(assimilation_efficiency_matrix)
        defined_parameters["P"]["assimilation_efficiency"] = Dict(
            k => interaction_args["P"][k] for
            k in ["can_eat", "can_be_eaten", "assimilation_efficiency"]
        )
        defined_parameters["Z"]["assimilation_efficiency"] = Dict(
            k => interaction_args["Z"][k] for
            k in ["can_eat", "can_be_eaten", "assimilation_efficiency"]
        )
    end

    emergent_parameters = compute_allometric_parameters(defined_parameters)

    if !isnothing(palatability_matrix)
        if !(size(palatability_matrix) == (n_phyto + n_zoo, n_phyto + n_zoo))
            throw(
                ArgumentError(
                    "palatability_matrix must have size $((n_phyto+n_zoo, n_phyto+n_zoo))"
                ),
            )
        end
        emergent_parameters["palatability_matrix"] = palatability_matrix
    end

    if !isnothing(assimilation_efficiency_matrix)
        if !(size(assimilation_efficiency_matrix) == (n_phyto + n_zoo, n_phyto + n_zoo))
            throw(
                ArgumentError(
                    "assimilation_efficiency_matrix must have size $((n_phyto+n_zoo, n_phyto+n_zoo))",
                ),
            )
        end
        emergent_parameters["assimilation_efficiency_matrix"] =
            assimilation_efficiency_matrix
    end

    # combine emergent parameters with remaining user defined parameters
    parameters = NamedTuple(
        Symbol(k) => v for (k, v) in merge(bgc_args, emergent_parameters)
    )

    # create tracer functions
    plankton_array = vcat(
        [Symbol("P$i") for i in 1:n_phyto], [Symbol("Z$i") for i in 1:n_zoo]
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
    return define_tracer_functions(parameters, tracers)
end

end # module
