# Lightweight diagnostic helpers.
module Diagnostics

import Oceananigans: time_step!

export box_model_budget
export box_model_mass_balance

include("box_model.jl")

end # module
