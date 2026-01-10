using Agate
using Test

include("test_biogeochemistry.jl")
include("test_constructor.jl")
include("test_parameters.jl")
include("test_library.jl")

# integration tests
include("test_box_model.jl")

# model checks
include("test_mass_balance.jl")
