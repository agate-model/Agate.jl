"""Runtime parameterization and box-ODE utilities."""
module Runtime

using Oceananigans.Biogeochemistry:
    AbstractContinuousFormBiogeochemistry,
    required_biogeochemical_auxiliary_fields,
    required_biogeochemical_tracers

export parameterized
export ode_problem
export active_parameters

include("active_parameters.jl")
include("ode_problem.jl")

end # module
