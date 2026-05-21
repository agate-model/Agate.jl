# Lightweight diagnostic helpers.
module Diagnostics

import Oceananigans: time_step!

export box_model_budget
export box_model_mass_balance
export ode_problem
export ode_initial_state

include("box_model.jl")
include("ode_problem.jl")

end # module
