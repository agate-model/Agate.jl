"""Scientific process authoring and setup-time normalization."""
module Processes

using ..Configuration: Population, Pool
using ..Factories:
    AbstractBGCFactory, default_components, default_processes, parameter_definitions

export AbstractProcess, AbstractFormulation
export Smith, Monod, PreferentialGrazing
export LinearMortality, QuadraticMortality, LinearRemineralization, PartitionRouting
export Growth, NutrientResponse, Grazing, Mortality, ProductRouting, Remineralization
export ModelDefinition, NormalizedModelDefinition, NamedProcess
export ParameterRequirementIdentity
export process_id, process_kind, formulation, formulation_tag
export participants, drivers, rate_axes, driver_identities, normalize_model

include("authoring.jl")
include("normalization.jl")

end # module Processes
