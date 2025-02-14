using Agate
using Test

# modules
include("test_biogeochemistry.jl")
include("test_constructor.jl")
include("test_parameters.jl")

# library
include("test_library_predation.jl")
include("test_library_mortality.jl")

# integration tests
include("test_box_model.jl")

# model checks
include("test_mass_balance.jl")
