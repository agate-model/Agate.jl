"""Reusable tracer tendency builders."""
module Tendencies

using ..Equations: CompiledEquation
using ..Library.Mortality: linear_loss, quadratic_loss
using ..Library.Photosynthesis: smith_growth, geider_growth
using ..Library.Remineralization: linear_remineralization
using ..Runtime: tendency_inputs

include("config.jl")
include("coupling.jl")
include("helpers.jl")
include("reductions.jl")

using .Reductions:
    grazing_unassimilated_loss_sum,
    grazing_loss_sum,
    grazing_gain_sum,
    growth_sum,
    mortality_loss_sum

include("phytoplankton.jl")
include("zooplankton.jl")
include("inorganic.jl")
include("organic_matter.jl")
include("detritus.jl")

export TendencyConfig,
    NutrientCoupling,
    nutrient_coupling,
    phytoplankton_tendency,
    zooplankton_tendency,
    inorganic_tendency,
    organic_matter_tendency,
    detritus_tendency

end # module Tendencies
