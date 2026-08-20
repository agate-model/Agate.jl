# to build docs from terminal: julia --project=. -e 'using Pkg; Pkg.test()'
using Agate
using Test

include(joinpath("LegacyTendencies", "Tendencies.jl"))
include("test_helpers.jl")
include(joinpath("NPZD", "tracers.jl"))

include("test_tracer_functions.jl")
include("test_classrefs_and_tracer_accessors.jl")
include("test_library.jl")
include("test_components.jl")
include("test_processes.jl")
include("test_process_compilation.jl")
include("test_growth_compilation.jl")
include("test_grazing_compilation.jl")
include("test_nipizd_process_compilation.jl")
include("test_extensibility_compilation.jl")
include("test_direct_model_definition_construction.jl")
include("test_nipizd_canonical_process_path.jl")
include("test_legacy_architecture_cleanup.jl")
include("test_box_model.jl")

include("test_models_construct.jl")
include("test_forwarddiff.jl")
include("test_active_parameters.jl")
include("test_enzyme.jl")
include("test_parameter_directory.jl")
include("test_mass_balance.jl")

include("test_introspection.jl")
include("test_recipe_serialization.jl")

include("test_biogeochemistry.jl")

include("test_cycle01_contracts.jl")
