using Agate
using Test

include("test_helpers.jl")

include("test_tracer_functions.jl")
include("test_classrefs_and_tracer_accessors.jl")
include("test_library.jl")
include("test_box_model.jl")

include("test_models_construct.jl")
include("test_parameter_directory.jl")
include("test_mass_balance.jl")

include("test_introspection.jl")
include("test_variants.jl")

include("test_biogeochemistry.jl")
