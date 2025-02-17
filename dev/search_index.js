var documenterSearchIndex = {"docs":
[{"location":"api/#Constructors","page":"API","title":"Constructors","text":"","category":"section"},{"location":"api/#Size-structured-NPZD-model","page":"API","title":"Size-structured NPZD model","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Agate.Constructors.NiPiZD.construct","category":"page"},{"location":"api/#Agate.Constructors.NiPiZD.construct","page":"API","title":"Agate.Constructors.NiPiZD.construct","text":"construct(;\n    n_phyto=2,\n    n_zoo=2,\n    phyto_diameters=Dict(\n        \"min_diameter\" => 2, \"max_diameter\" => 10, \"splitting\" => \"log_splitting\"\n    ),\n    zoo_diameters=Dict(\n        \"min_diameter\" => 20, \"max_diameter\" => 100, \"splitting\" => \"linear_splitting\"\n    ),\n    phyto_args=DEFAULT_PHYTO_ARGS,\n    zoo_args=DEFAULT_ZOO_ARGS,\n    interaction_args=DEFAULT_INTERACTION_ARGS,\n    bgc_args=DEFAULT_BGC_ARGS,\n    palatability_matrix=nothing,\n    assimilation_efficiency_matrix=nothing,\n) -> DataType\n\nConstruct a size-structured NPZD model abstract type.\n\nThis constructor builds a size-structured plankton model with two plankton functional types: phytoplankton (P) and zooplankton (Z), each of which can be specified to have any number of size classes (n_phyto and n_zoo). In addition to plankton, the constructor implements idealized detritus (D) and nutrient (N) cycling by default, although more complex N and D cycling can also be defined using the nutrient_dynamics and detritus_dynamics arguments.\n\nDuring model construction, the size of each plankton determines photosynthetic growth rates, nutrient half saturation constants, predation rates, and predator-prey assimilation and palatability values. Alternatively, if manually defined predator-prey assimilation and palatability values are desired, these can be specified using the palatability_matrix and assimilation_efficiency_matrix arguments.\n\nNote that if non-default *_dynamics expressions are passed, the relevant *_args also need to be specified.\n\nThe type specification includes a photosynthetic active radiation (PAR) auxiliary field.\n\nKeywords\n\nn_phyto: number of phytoplankton in the model\nn_zoo: number of zooplankton in the model\nphyto_diameters: dictionary from which phyto diameters can be computed or a list of   values to use (as many as the model expects)\nzoo_diameters: dictionary from which zoo diameters can be computed or a list of   values to use (as many as the model expects)\nnutrient_dynamics: expression describing how nutrients change over time, see   Agate.Models.Tracers\ndetritus_dynamics: expression describing how detritus evolves over time, see   Agate.Models.Tracers\nphyto_dynamics: expression describing how phytoplankton grow, see Agate.Models.Tracers\nzoo_dynamics: expression describing how zooplankton grow, see Agate.Models.Tracers\nphyto_args: Dictionary of phytoplankton parameters, for default values see   Agate.Models.Constructors.DEFAULT_PHYTO_ARGS\nzoo_args: Dictionary of zooplankton parameters, for default values see   Agate.Models.Constructors.DEFAULT_ZOO_ARGS\ninteraction_args: Dictionary of arguments from which a palatability and assimilation  efficiency matrix between all plankton can be computed, for default values see   Agate.Models.Constructors.DEFAULT_INTERACTION_ARGS\nbgc_args: Dictionary of biogeochemistry parameters related to nutrient and detritus, for   default values see Agate.Models.Constructors.DEFAULT_BGC_ARGS\npalatability_matrix: optional palatability matrix passed as a NamedArray, if provided   then interaction_args are not used to compute this\nassimilation_efficiency_matrix: optional assimilation efficiency matrix passed as a   NamedArray, if provided then interaction_args are not used to compute this\n\nExample\n\nusing Agate.Constructors: NiPiZD\n\nn2p2zd = NiPiZD.construct()\nn2p2zd_model_obj = n2p2zd()\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"API","title":"API","text":"Agate.Constructors.NiPiZD.instantiate","category":"page"},{"location":"api/#Agate.Constructors.NiPiZD.instantiate","page":"API","title":"Agate.Constructors.NiPiZD.instantiate","text":"instantiate(\n    bgc_type;\n    phyto_diameters=Dict(\n        \"min_diameter\" => 2, \"max_diameter\" => 10, \"splitting\" => \"log_splitting\"\n    ),\n    zoo_diameters=Dict(\n        \"min_diameter\" => 20, \"max_diameter\" => 100, \"splitting\" => \"linear_splitting\"\n    ),\n    phyto_args=DEFAULT_PHYTO_ARGS,\n    zoo_args=DEFAULT_ZOO_ARGS,\n    interaction_args=DEFAULT_INTERACTION_ARGS,\n    bgc_args=DEFAULT_BGC_ARGS,\n    palatability_matrix=nothing,\n    assimilation_efficiency_matrix=nothing,\n)\n\nA function to instantiate an object of NiPiZD.construct() model type.\n\nThe type specifies the number of phytoplankton and zooplankton in the model and includes default parameter values. The instantiate method can be used to override the default values of any of the model parameters or plankton diameters.\n\nArguments\n\nbgc_type: subtype of Oceananigans.Biogeochemistry returned by NiPiZD.construct()  with a specified number of phytoplankton and zooplankton\n\nKeywords\n\nphyto_diameters: dictionary from which phyto diameters can be computed or a list of\n\nvalues to use (as many as the model expects)\n\nzoo_diameters: dictionary from which zoo diameters can be computed or a list of   values to use (as many as the model expects)\nnutrient_dynamics: expression describing how nutrients change over time, see   Agate.Models.Tracers\ndetritus_dynamics: expression describing how detritus evolves over time, see   Agate.Models.Tracers\nphyto_dynamics: expression describing how phytoplankton grow, see Agate.Models.Tracers\nzoo_dynamics: expression describing how zooplankton grow, see Agate.Models.Tracers\nphyto_args: Dictionary of phytoplankton parameters, for default values see   Agate.Models.Constructors.DEFAULT_PHYTO_ARGS\nzoo_args: Dictionary of zooplankton parameters, for default values see   Agate.Models.Constructors.DEFAULT_ZOO_ARGS\ninteraction_args: Dictionary of arguments from which a palatability and assimilation  efficiency matrix between all plankton can be computed, for default values see   Agate.Models.Constructors.DEFAULT_INTERACTION_ARGS\nbgc_args: Dictionary of constant parameters used in growth functions (i.e., not size   dependant plankton parameters as well as biogeochemistry parameters related to nutrient   and detritus, for default values see Agate.Models.Constructors.DEFAULT_CONSTANT_ARGS\npalatability_matrix: optional palatability matrix passed as a NamedArray, if provided   then interaction_args are not used to compute this\nassimilation_efficiency_matrix: optional assimilation efficiency matrix passed as a   NamedArray, if provided then interaction_args are not used to compute this\n\nExample\n\nusing Agate.Constructors: NiPiZD\n\nn2p2zd = NiPiZD.construct()\n\n# change some parameter values\nphyto_args = NiPiZD.DEFAULT_PHYTO_ARGS\nphyto_args[\"allometry\"][\"maximum_growth_rate\"][\"a\"] = 2\nn2p2zd_model_obj = NiPiZD.instantiate(n2p2zd; phyto_args=phyto_args)\n\n\n\n\n\n","category":"function"},{"location":"api/#Low-level-API","page":"API","title":"Low level API","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Agate.Models.Biogeochemistry.define_tracer_functions","category":"page"},{"location":"api/#Agate.Models.Biogeochemistry.define_tracer_functions","page":"API","title":"Agate.Models.Biogeochemistry.define_tracer_functions","text":"define_tracer_functions(\n    parameters,\n    tracers;\n    auxiliary_fields=[:PAR],\n    helper_functions=nothing,\n    sinking_tracers=nothing,\n    grid=nothing,\n    open_bottom=false,\n) -> DataType\n\nCreate an Oceananigans.Biogeochemistry model type.\n\nArguments\n\nparameters: named sequence of values of the form ((<field name> = <default value>, ...)\ntracers: dictionary of the form (<name> => <expression>, ...)\n\nKeywords\n\nauxiliary_fields: an iterable of auxiliary field variables, defaults to [:PAR,]\nhelper_functions: optional path to a file of helper functions used in tracer expressions\nsinking_tracers: optional NamedTuple of sinking speeds (passed as positive values) of  the form (<tracer name> = <speed>, ...)\ngrid: optional Oceananigans grid object defining the geometry to build the model on, must  be passed if sinking_tracers is defined\nopen_bottom: indicates whether the sinking velocity should be smoothly brought to zero  at the bottom to prevent the tracers leaving the domain, defaults to true, which means  the bottom is open and the tracers leave (i.e., no slowing of velocity to 0 is applied)\n\nNote that the field names defined in parameters can't be any of [:x, :y, :z, :t], as these are reserved for coordinates, and they must include all parameters used in the tracers expressions. The expressions must use methods that are either defined within this module or passed in the helper_functions file.\n\nExample\n\nusing Agate\n\nparameters = (α=2 / 3, β=4 / 3, δ=1, γ=1)\ntracers = Dict(\"R\" => :(α * R - β * R * F), \"F\" => :(-γ * F + δ * R * F))\nLV = define_tracer_functions(parameters, tracers)\n\n\n\n\n\n","category":"function"},{"location":"library/#Allometry","page":"Library","title":"Allometry","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"Modules=[Agate.Library.Allometry]","category":"page"},{"location":"library/#Agate.Library.Allometry.allometric_palatability_unimodal-Tuple{Dict, Dict}","page":"Library","title":"Agate.Library.Allometry.allometric_palatability_unimodal","text":"allometric_palatability_unimodal(prey_data, predator_data)\n\nCalculates the unimodal allometric palatability of prey based on predator-prey diameters.\n\nThis function extracts prey_diameter, predator_diameter, optimum_predator_prey_ratio,  and specificity from the provided dictionaries and calculates the palatability using the diameter-based formula.\n\nNote that this formulation differs from the currently operational MITgcm-DARWIN model as it uses diameter instead of volumes and is structurally different. However, both formulations result in a unimodal response the width and optima are modulated by the optimumpredatorprey ratio and the specificity.\n\nArguments\n\nprey_data: A dictionary containing prey-specific data:\ndiameters: Diameter of the prey.\npredator_data: A dictionary containing predator-specific data:\ncan_eat: A binary value (1 or 0) indicating if the predator can consume prey. If this is set to 0, palatability is set to 0.\ndiameters: Diameter of the predator.\noptimum_predator_prey_ratio: The optimal predator-prey diameter ratio for the predator.\nspecificity: A parameter controlling how sharply the palatability decreases away from the optimal ratio.\n\nReturns\n\npalatability: A number between 0 and 1 representing the palatability.\n\n\n\n\n\n","category":"method"},{"location":"library/#Agate.Library.Allometry.allometric_palatability_unimodal_protection-Tuple{Dict, Dict}","page":"Library","title":"Agate.Library.Allometry.allometric_palatability_unimodal_protection","text":"allometric_palatability_unimodal_protection(prey_data, predator_data)\n\nCalculates the unimodal allometric palatability of prey, accounting for additional prey protection mechanisms.\n\nThe function uses a modified unimodal relationship defined by: palatability = prey_protection / (1 + (predator_prey_ratio - predator_prey_optimum)^2)^predator_specificity\n\nArguments\n\nprey_data: A dictionary containing prey-specific data:\ndiameters: Diameter of the prey.\nprotection: A scaling factor between 0 and 1 representing additional protection mechanisms of the prey.\npredator_data: A dictionary containing predator-specific data:\ncan_eat: A binary value (1 or 0) indicating if the predator can consume prey. If this is set to 0, palatability is set to 0.\ndiameters: Diameter of the predator.\noptimum_predator_prey_ratio: The optimal predator-prey diameter ratio for the predator.\nspecificity: A parameter controlling how sharply the palatability decreases away from the optimal ratio.\n\nReturns\n\npalatability: A number between 0 and prey_protection representing the palatability.\n\n\n\n\n\n","category":"method"},{"location":"library/#Agate.Library.Allometry.allometric_scaling_power-Tuple{Number, Number, Number}","page":"Library","title":"Agate.Library.Allometry.allometric_scaling_power","text":"allometric_scaling_power(a, b, d::Number)\n\nAllometric scaling function using the power law for cell volume.\n\nArguments\n\na: scale\nb: exponent\nd: cell equivalent spherical diameter (ESD)\n\n\n\n\n\n","category":"method"},{"location":"library/#Growth","page":"Library","title":"Growth","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"Modules=[Agate.Library.Growth]","category":"page"},{"location":"library/#Agate.Library.Growth.default_PC-NTuple{5, Any}","page":"Library","title":"Agate.Library.Growth.default_PC","text":"PC = PCᵐᵃˣ * γⁿᵘᵗ *  γˡⁱᵍʰᵗ * fᵗᵉᵐᵖ *  γᶜᵒ²\n\nCarbon-specific growth rate for plankton (Default MITgcm-DARWIN formulation).\n\nArguments\n\nPCᵐᵃˣ: maximum carbon-specific growth rate\nγⁿᵘᵗ: nutrient limition\nγˡⁱᵍʰᵗ: light limition\nfᵗᵉᵐᵖ: temperature limitation\nγᶜᵒ²: carbon dioxide limitation\n\n\n\n\n\n","category":"method"},{"location":"library/#Mortality","page":"Library","title":"Mortality","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"Modules=[Agate.Library.Mortality]","category":"page"},{"location":"library/#Agate.Library.Mortality.linear_loss-Tuple{Any, Any}","page":"Library","title":"Agate.Library.Mortality.linear_loss","text":"Linear mortality rate. In this formulation mortality is constant, and can be interpreted as a \"closure term\" for low density predation and and other death terms.\n\nArguments\n\nP: plankton concentration\nl: mortality rate\n\n\n\n\n\n","category":"method"},{"location":"library/#Agate.Library.Mortality.net_linear_loss-Tuple{Any, Any, Any}","page":"Library","title":"Agate.Library.Mortality.net_linear_loss","text":"Net loss of all plankton due to linear mortality.\n\nArguments\n\nP: NamedArray which includes all plankton concentration values\nlinear_mortality: NamedArray of plankton linear mortality rates\n\n\n\n\n\n","category":"method"},{"location":"library/#Agate.Library.Mortality.net_quadratic_loss","page":"Library","title":"Agate.Library.Mortality.net_quadratic_loss","text":"Net loss of all plankton due to quadratic mortality.\n\nArguments\n\nP: NamedArray which includes all plankton concentration values\nquadratic_mortality: plankton quadratic mortality rate\nplankton_type_prefix: Array of prefixes used in plankton names to indicate their type,   use here to sum over only the relevant plankton (e.g., \"Z\" for zooplankton)\n\n\n\n\n\n","category":"function"},{"location":"library/#Agate.Library.Mortality.quadratic_loss-Tuple{Any, Any}","page":"Library","title":"Agate.Library.Mortality.quadratic_loss","text":"Quadratic mortality coefficient. In this formulation mortality increases exponentially with plankton biomass and is often interpreted to represent viral processes and non-represented density-dependent predation effects.\n\nArguments\n\nP: plankton concentration\nl: mortality rate\n\n\n\n\n\n","category":"method"},{"location":"library/#Nutrients","page":"Library","title":"Nutrients","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"Modules=[Agate.Library.Nutrients]","category":"page"},{"location":"library/#Agate.Library.Nutrients.monod_limitation-Tuple{Any, Any}","page":"Library","title":"Agate.Library.Nutrients.monod_limitation","text":"R / (kᵣ + R)\n\nMonod formulation of nutrient limitation, which is based on Michaelis-Menten enzyme kinetics.\n\nArguments\n\nR: nutrient (e.g. N, P, Si)\nkᵣ: nutrient half saturation constant\n\nNote that sometimes this formulation is also used for Predation ('Holling type 2').\n\n\n\n\n\n","category":"method"},{"location":"library/#Photosynthesis","page":"Library","title":"Photosynthesis","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"Modules=[Agate.Library.Photosynthesis]","category":"page"},{"location":"library/#Agate.Library.Photosynthesis.light_limitation_geider-NTuple{4, Any}","page":"Library","title":"Agate.Library.Photosynthesis.light_limitation_geider","text":"Pᶜₘₐₓ[1-exp((-αᶜʰˡθᶜE₀)/Pᶜₘₐₓ)]\n\nA light limitation function which depends on the cellular ratio of chlorophyll to carbon. This formulation is based on equation (4) from Geider et al., 1998.\n\nArguments\n\nPAR: photosynthetic active radiation (E₀)\nmaximum_growth_rate: maximum growth rate before nutrient limitation (Pᶜₘₐₓ)\nphotosynthetic_slope: initial photosynthetic slope (αᶜʰˡ)\nchlorophyll_to_carbon_ratio: ratio between cellular chlorophyll and carbon (θᶜ)\n\n\n\n\n\n","category":"method"},{"location":"library/#Agate.Library.Photosynthesis.light_limitation_smith-Tuple{Any, Any, Any}","page":"Library","title":"Agate.Library.Photosynthesis.light_limitation_smith","text":"α * PAR / sqrt(μ₀ ^ 2 + α ^ 2 * PAR ^ 2)\n\nSmith 1936 formulation of light limitation (also see Evans and Parslow, 1985).\n\nArguments\n\nPAR: photosynthetic active radiation\nα: initial photosynthetic slope\nμ₀: maximum growth rate at T = 0 °C (this seems weird?, from Kuhn 2015)\n\n\n\n\n\n","category":"method"},{"location":"library/#Agate.Library.Photosynthesis.net_photosynthetic_growth_single_nutrient","page":"Library","title":"Agate.Library.Photosynthesis.net_photosynthetic_growth_single_nutrient","text":"Net photosynthetic growth of all plankton assuming geider light limitation.\n\nArguments\n\nN: Nutrient\nP: NamedArray which includes all plankton concentration values\nPAR: PAR\nmaximum_growth_rate: NamedArray of all plankton maximum growth rates\nnutrient_half_saturation: NamedArray of all plankton nutrient half saturation constants\nalpha: initial photosynthetic slope\nplankton_type_prefix: Array of prefixes used in plankton names to indicate their type,   use here to sum over only the relevant plankton (e.g., \"P\" for phytoplankton)\n\n\n\n\n\n","category":"function"},{"location":"library/#Agate.Library.Photosynthesis.net_photosynthetic_growth_single_nutrient_geider_light","page":"Library","title":"Agate.Library.Photosynthesis.net_photosynthetic_growth_single_nutrient_geider_light","text":"Net photosynthetic growth of all plankton.\n\nArguments\n\nN: Nutrient\nP: NamedArray which includes all plankton concentration values\nPAR: PAR\nmaximum_growth_rate: NamedArray of all plankton maximum growth rates\nnutrient_half_saturation: NamedArray of all plankton nutrient half saturation constants\nphotosynthetic_slope: initial photosynthetic slope (αᶜʰˡ)\nchlorophyll_to_carbon_ratio: ratio between cellular chlorophyll and carbon (θᶜ)\nplankton_type_prefix: Array of prefixes used in plankton names to indicate their type,   use here to sum over only the relevant plankton (e.g., \"P\" for phytoplankton)\n\n\n\n\n\n","category":"function"},{"location":"library/#Agate.Library.Photosynthesis.photosynthetic_growth_single_nutrient-NTuple{6, Any}","page":"Library","title":"Agate.Library.Photosynthesis.photosynthetic_growth_single_nutrient","text":"Single nutrient monod smith photosynthetic growth (used, for example, in Kuhn 2015).\n\nArguments\n\nN: nutrient concentration\nP: phytoplankton concentration\nPAR: photosynthetic active radiation\nμ₀: maximum growth rate at T = 0 °C\nkₙ: nutrient half saturation\nα: initial photosynthetic slope\n\n\n\n\n\n","category":"method"},{"location":"library/#Agate.Library.Photosynthesis.photosynthetic_growth_single_nutrient_geider_light-NTuple{7, Any}","page":"Library","title":"Agate.Library.Photosynthesis.photosynthetic_growth_single_nutrient_geider_light","text":"Single nutrient geider photosynthetic growth.\n\nArguments\n\nN: nutrient concentration\nP: phytoplankton concentration\nPAR: photosynthetic active radiation\nmaximum_growth_rate: maximum growth rate before nutrient limitation (Pᶜₘₐₓ)\nkₙ: nutrient half saturation\nphotosynthetic_slope: initial photosynthetic slope (αᶜʰˡ)\nchlorophyll_to_carbon_ratio: ratio between cellular chlorophyll and carbon (θᶜ)\n\n\n\n\n\n","category":"method"},{"location":"library/#Agate.Library.Photosynthesis.γˡⁱᵍʰᵗ-NTuple{4, Any}","page":"Library","title":"Agate.Library.Photosynthesis.γˡⁱᵍʰᵗ","text":"γˡⁱᵍʰᵗ = (1 - ℯ^(kˢᵃᵗ*I)) * ℯ^kⁱⁿʰ * nˡⁱᵍʰᵗ\n\nLight limitation for plankton (Default MITgcm-DARWIN formulation).\n\nArguments\n\nI: irradiance\nkˢᵃᵗ:  half saturation constant of light saturation\nkⁱⁿʰ: half saturation constant of light inhibition\nnˡⁱᵍʰᵗ: light penalty term\n\n\n\n\n\n","category":"method"},{"location":"library/#Predation","page":"Library","title":"Predation","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"Modules=[Agate.Library.Predation]","category":"page"},{"location":"library/#Agate.Library.Predation.assimilation_efficiency_emergent_binary-Tuple{Any, Any}","page":"Library","title":"Agate.Library.Predation.assimilation_efficiency_emergent_binary","text":"assimilation_efficiency_emergent_binary(prey_data, predator_data)\n\nDetermines the assimilation efficiency of a predator consuming prey, based on binary conditions of edibility.\n\nThe function evaluates whether the predator can eat the prey and whether the prey can be consumed, and assigns the assimilation efficiency accordingly.\n\nArguments\n\nprey_data: A dictionary containing prey-specific data:\ncan_be_eaten: A binary value (1 or 0) indicating if the prey can be consumed by the predator.\npredator_data: A dictionary containing predator-specific data:\ncan_eat: A binary value (1 or 0) indicating if the predator can consume prey.\nassimilation_efficiency: The efficiency with which the predator assimilates nutrients from the prey if the conditions are met.\n\nReturns\n\nassimilation_efficiency:\nIf can_eat is 1 and can_be_eaten is 1, returns the predator's assimilation_efficiency.\nOtherwise, returns 0.\n\n\n\n\n\n","category":"method"},{"location":"library/#Agate.Library.Predation.holling_type_2-Tuple{Real, Real}","page":"Library","title":"Agate.Library.Predation.holling_type_2","text":"Holling's \"type II\" functional response as describe in Holling 1959. The function is similar to the Monod equation and Michaelis-Menten equation of for enzyme kinetics. The formulation is characterized by decelerating predation as prey concentrations increase.\n\nArguments\n\nR: prey density\nk: prey density at which predation is half it's maximum rate\n\n\n\n\n\n","category":"method"},{"location":"library/#Agate.Library.Predation.net_predation_assimilation_loss_preferential","page":"Library","title":"Agate.Library.Predation.net_predation_assimilation_loss_preferential","text":"Net predator assimilation loss of all plankton.\n\nArguments\n\nP: NamedArray which includes all plankton concentration values\nholling_half_saturation: NamedArray of all plankton predation half saturation constants\nmaximum_predation_rate: NamedArray of all plankton maximum predation rates\npalatability: NamedArray of all plankton palatabilities where:\neach row is a predator\neach column is a prey\nvalues are accessed as palat[predator, prey]\nfor a non-predator [i,:]=0\nassimilation_efficiency: NamedArray of all plankton assimilation efficiencies where:\neach row is a predator\neach column is a prey\nvalues are accessed as palat[predator, prey]\nfor a non-predator [i,:]=0\nplankton_type_prefix: Array of prefixes used in plankton names to indicate their type,   use here to sum over only predator plankton (e.g., \"Z\" for zooplankton)\n\n\n\n\n\n","category":"function"},{"location":"library/#Agate.Library.Predation.predation_assimilation_loss_idealized-NTuple{5, Any}","page":"Library","title":"Agate.Library.Predation.predation_assimilation_loss_idealized","text":"Estimates the rate at which plankton predation gain is lost to the environment due to inefficient assimilation efficiency (e.g. 'sloppy feeding').\n\nNote that this differs from the predationgainidealized as this function represents the transfer of biomass from the prey to the environment rather than the transfer of biomass from the prey to the predator.\n\nArguments\n\nP: phytoplankton concentration\nZ: zooplankton concentration\nβ: assimilation efficiency of prey to the predator\ngₘₐₓ: maximum grazing rate of the predator\nkₚ: grazing/holling half saturation\n\n\n\n\n\n","category":"method"},{"location":"library/#Agate.Library.Predation.predation_assimilation_loss_preferential-NTuple{6, Any}","page":"Library","title":"Agate.Library.Predation.predation_assimilation_loss_preferential","text":"Estimates the rate at which plankton predation gain is lost to the environment due to inefficient assimilation efficiency (e.g. 'sloppy feeding').\n\nNote that this differs from the predationgainpreferential as this function represents the transfer of biomass from the prey to the environment rather than the transfer of biomass from the prey to the predator.\n\nArguments\n\nP: phytoplankton concentration\nZ: zooplankton concentration\nβ: assimilation efficiency of prey to the predator\ngₘₐₓ: maximum grazing rate of the predator\nkₚ: grazing/holling half saturation\npalatability: the likelihood at which the predator feeds on the prey\n\n\n\n\n\n","category":"method"},{"location":"library/#Agate.Library.Predation.predation_gain_idealized-NTuple{5, Any}","page":"Library","title":"Agate.Library.Predation.predation_gain_idealized","text":"Estimates the gain rate of Z (predator) feeding on P (prey). In this formulation predation rate is multiplied by an assimilation efficiency (β) which represents 'sloppy feeding'.\n\nArguments\n\nP: phytoplankton concentration\nZ: zooplankton concentration\nβ: assimilation efficiency\ngₘₐₓ: maximum grazing rate\nkₚ: grazing/holling half saturation\n\n\n\n\n\n","category":"method"},{"location":"library/#Agate.Library.Predation.predation_gain_preferential-NTuple{6, Any}","page":"Library","title":"Agate.Library.Predation.predation_gain_preferential","text":"Estimates the gain rate of Z (predator) feeding on P (prey). In this formulation predation rate is multiplied by an assimilation efficiency (β) which represents 'sloppy feeding'.\n\nArguments\n\nP: phytoplankton concentration\nZ: zooplankton concentration\nβ: assimilation efficiency\ngₘₐₓ: maximum grazing rate\nkₚ: grazing/holling half saturation\npalatability: the likelihood at which the predator feeds on the prey\n\n\n\n\n\n","category":"method"},{"location":"library/#Agate.Library.Predation.predation_loss_idealized-NTuple{4, Any}","page":"Library","title":"Agate.Library.Predation.predation_loss_idealized","text":"Estimates the loss rate of P (prey), to Z (predator). In this formulation predator-prey interactions are modulated both by their density (Holling type 2) and the prey-predator palatability.\n\nArguments\n\nP: phytoplankton concentration\nZ: zooplankton concentration\ngₘₐₓ: maximum grazing rate\nkₚ: prey density at which predation is half it's maximum rate\n\n\n\n\n\n","category":"method"},{"location":"library/#Agate.Library.Predation.predation_loss_preferential-NTuple{5, Any}","page":"Library","title":"Agate.Library.Predation.predation_loss_preferential","text":"Estimates the loss rate of P (prey), to Z (predator). In this formulation predator-prey interactions are modulated both by their density (Holling type 2) and the prey-predator palatability.\n\nArguments\n\nP: phytoplankton concentration\nZ: zooplankton concentration\ngₘₐₓ: maximum grazing rate\nkₚ: prey density at which predation is half it's maximum rate\npalatability: the likelihood at which the predator feeds on the prey\n\n\n\n\n\n","category":"method"},{"location":"library/#Agate.Library.Predation.summed_predation_assimilation_loss_preferential-NTuple{6, Any}","page":"Library","title":"Agate.Library.Predation.summed_predation_assimilation_loss_preferential","text":"Estimates the total assimilation loss of the predator (P[predator_name]) feeding on all plankton.\n\nFor plankton P[predator_name], the function loops over each prey (P[prey_name]) to estimate the total assimilation loss during predation.\n\nArguments\n\npredator_name: name of the predator, e.g. P[predator_name]\nP: NamedArray which includes all plankton concentration values\nmaximum_predation_rate: NamedArray of all plankton predation rates\nholling_half_saturation: plankton predation half saturation constant\npalatability: NamedArray of all plankton palatabilities where:\neach row is a predator\neach column is a prey\nvalues are accessed as palat[predator, prey]\nfor a non-predator [i,:]=0\nassimilation_efficiency: NamedArray of all plankton assimilation efficiencies where:\neach row is a predator\neach column is a prey\nvalues are accessed as palat[predator, prey]\nfor a non-predator [i,:]=0\n\n\n\n\n\n","category":"method"},{"location":"library/#Agate.Library.Predation.summed_predation_gain_preferential-NTuple{6, Any}","page":"Library","title":"Agate.Library.Predation.summed_predation_gain_preferential","text":"Estimates the total predation gain of the predator (P[predator_name]) feeding on all plankton.\n\nFor plankton P[predator_name], the function loops over each prey (P[prey_name]) to estimate the total gain due to predation.\n\nArguments\n\npredator_name: name of the predator, e.g. P[predator_name]\nP: NamedArray which includes all plankton concentration values\nmaximum_predation_rate: NamedArray of all plankton predation rates\nholling_half_saturation: predation half saturation constant\npalatability: NamedArray of all plankton palatabilities where:\neach row is a predator\neach column is a prey\nvalues are accessed as palat[predator, prey]\nfor a non-predator [i,:]=0\nassimilation_efficiency: NamedArray of all plankton assimilation efficiencies where:\neach row is a predator\neach column is a prey\nvalues are accessed as palat[predator, prey]\nfor a non-predator [i,:]=0\n\n\n\n\n\n","category":"method"},{"location":"library/#Agate.Library.Predation.summed_predation_loss_preferential","page":"Library","title":"Agate.Library.Predation.summed_predation_loss_preferential","text":"Estimates the total loss rate of the prey P[prey_name] to predation.\n\nFor plankton P[prey_name], the function loops over each predator to estimate the total loss of plankton prey_name due to predation.\n\nArguments\n\nprey_name: name of the prey plankton to access value as P[prey_name]\nP: NamedArray which includes all plankton\nmaximum_predation_rate: NamedArray of all plankton predation rates\nholling_half_saturation: predation half saturation constant\npalatability: NamedArray of all plankton palatabilities where:\neach row is a predator\neach column is a prey\nvalues are accessed as palat[predator, prey]\nfor a non-predator [i,:]=0\nplankton_type_prefix: Array of prefixes used in plankton names to indicate their type,   use here to sum over only predator plankton (e.g., \"Z\" for zooplankton)\n\n\n\n\n\n","category":"function"},{"location":"library/#Remineralization","page":"Library","title":"Remineralization","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"Modules=[Agate.Library.Remineralization]","category":"page"},{"location":"library/#Agate.Library.Remineralization.remineralization_idealized-Tuple{Any, Any}","page":"Library","title":"Agate.Library.Remineralization.remineralization_idealized","text":"Idealized remineralization of detritus into dissolved nutrients.\n\nArguments\n\nD: detritus\nr: remineralization rate\n\n\n\n\n\n","category":"method"},{"location":"#Agate.jl","page":"About","title":"Agate.jl","text":"","category":"section"},{"location":"","page":"About","title":"About","text":"A Julia library to build flexible and composable aquatic ecosystems.","category":"page"}]
}
