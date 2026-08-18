# # [Defining your own model] (@id external_model_family)

# Agate provides the [Agate.jl-NiPiZD](@ref NiPiZD) model, but it is also designed so that you can define your
# own model families. Here we define a simple one-tracer model in which a tracer decays
# at a constant rate, using the same factory interface available to external packages.

using Agate.Construction: construct_factory
using Agate.Equations: CompiledEquation
using Agate.Factories: AbstractBGCFactory, ParameterDefinition, ParameterSpec, ConstDefault
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers
using Oceananigans.Units: day

import Agate.Construction: recipe_family, recipe_factory
import Agate.Factories:
    parameter_definitions,
    default_community,
    default_plankton_dynamics,
    default_biogeochem_dynamics

nothing #hide

# ## Define the family

struct DecayFactory <: AbstractBGCFactory end

# These methods give the family a stable recipe identifier and let Agate reconstruct its
# factory when a recipe is loaded.

recipe_family(::DecayFactory) = :Decay
recipe_factory(::Val{:Decay}) = DecayFactory()
nothing #hide

# The model has one scalar parameter, with a default decay rate of 0.1 day⁻¹.

parameter_definitions(::DecayFactory) = (
    ParameterDefinition(ParameterSpec(:decay_rate, :scalar), ConstDefault(0.1 / day)),
)
nothing #hide

# There are no plankton groups in this model.

default_community(::DecayFactory) = (;)
default_plankton_dynamics(::DecayFactory) = (;)
nothing #hide

# The single tracer follows ``dX/dt = -kX``.

decay_tendency() = CompiledEquation(
    (bgc, x, y, z, t, X) -> -bgc.parameters.decay_rate * X
)

default_biogeochem_dynamics(::DecayFactory) = (X=decay_tendency,)

nothing #hide

# ## Construct the model

# `X` does not require any auxiliary fields, so we disable Agate's default `PAR` field.

bgc = construct_factory(DecayFactory(); auxiliary_fields=())

println(required_biogeochemical_tracers(bgc))
println("decay rate = ", bgc.parameters.decay_rate * day, " day⁻¹")
println("tendency at X = 2: ", bgc(Val(:X), 0, 0, 0, 0, 2.0) * day, " day⁻¹")

nothing #hide
