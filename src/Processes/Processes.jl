"""Scientific process authoring and setup-time normalization."""
module Processes

using ..Configuration: Population, Pool
using ..Factories:
    AbstractBGCFactory, default_components, default_processes, parameter_definitions
using ..Library.Mortality: linear_loss, quadratic_loss

export AbstractProcess, AbstractFormulation
export Smith, Monod, PreferentialGrazing
export LinearMortality, QuadraticMortality, LinearRemineralization, PartitionRouting
export Growth, NutrientResponse, Grazing, Mortality, ProductRouting, Remineralization
export ModelDefinition, NormalizedModelDefinition, NamedProcess
export ParameterRequirementIdentity
export process_id, process_kind, formulation, formulation_tag
export participants, drivers, rate_axes, driver_identities, normalize_model
export process_rate

include("authoring.jl")
include("rates.jl")
include("normalization.jl")

end # module Processes
