# # [External model family] (@id external_model_family_example)

# Agate provides the [Agate.jl-NiPiZD](@ref NiPiZD) model, but it is also designed so that you can define your
# own model families. Here we define a simple one-tracer model in which a tracer decays
# at a constant rate, using the same factory interface available to external packages.

# ## Loading dependencies

using Agate.Construction: construct_factory
using Agate.Equations: CompiledEquation
using Agate.Factories:
    AbstractBGCFactory, ParameterDefinition, ParameterSpec, ConstDefault
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers
using Oceananigans.Units: day

import Agate.Construction: recipe_family, recipe_factory
import Agate.Factories:
    parameter_definitions,
    default_community,
    default_plankton_dynamics,
    default_biogeochem_dynamics

nothing #hide

# ## Define the model family

# Agate represents each model family with a *factory*. A Julia `struct` defines a new type.
# `DecayFactory` has no fields because there is no information to store in the factory itself;
# its type is enough to tell Agate which model family we are working with. Subtyping
# `AbstractBGCFactory` marks it as an Agate biogeochemistry factory.
#
# Using a type in this way lets Agate use Julia's *multiple dispatch*. For example,
# `parameter_definitions(::DecayFactory)` below means "use this method of
# `parameter_definitions` whenever the argument is a `DecayFactory`". A more complicated
# factory could also have fields containing information needed to define the family.

struct DecayFactory <: AbstractBGCFactory end

# A `Symbol` is a lightweight Julia identifier written with a leading colon, such as
# `:Decay`, `:decay_rate`, or `:X`. Agate uses symbols for stable programmatic names: here
# `:Decay` is the name stored in a recipe to identify this model family.
#
# `Val{:Decay}` moves that symbol into the *type* of an argument. This lets Julia dispatch
# directly on the identifier. When Agate reads a recipe whose family is `:Decay`, it can call
# `recipe_factory(Val(:Decay))`, and Julia selects the method below. This avoids a central
# list or `if` statement that would need to be edited whenever an external package adds a
# model family.

recipe_family(::DecayFactory) = :Decay
recipe_factory(::Val{:Decay}) = DecayFactory()

nothing #hide

# ## Define the model parameter

# Parameters are declared independently of the tracer equations. Here the model has one
# scalar decay rate with a default value of 0.1 day⁻¹.
#
# `ParameterSpec(:decay_rate, :scalar)` says that the parameter is named `:decay_rate` and
# contains one scalar value. `ConstDefault` supplies its default. The method is written for
# `::DecayFactory`, so these parameter definitions belong specifically to the Decay family.

function parameter_definitions(::DecayFactory)
    return (
        ParameterDefinition(
            ParameterSpec(:decay_rate, :scalar),
            ConstDefault(0.1 / day),
        ),
    )
end

nothing #hide

# ## Define the tracer equation

# This family has no plankton groups, so its community and plankton dynamics are empty.
# The empty named tuple `(; )` tells Agate that there is nothing to configure for those
# parts of this model.

default_community(::DecayFactory) = (;)
default_plankton_dynamics(::DecayFactory) = (;)

# The model has one tracer, `X`, which decays according to
#
# ```math
# \frac{dX}{dt} = -kX,
# ```
#
# where `k` is `decay_rate`. `CompiledEquation` wraps the Julia function that evaluates this
# tendency. Agate passes the model (`bgc`), position (`x`, `y`, `z`), time (`t`), and tracer
# concentration (`X`) to the function. This simple equation only needs `bgc` and `X`.

decay_tendency() = CompiledEquation(
    (bgc, x, y, z, t, X) -> -bgc.parameters.decay_rate * X
)

# A named tuple connects tracer names to their tendency builders. `X=decay_tendency` means
# that the tracer identified by the symbol `:X` uses `decay_tendency` as its equation.

default_biogeochem_dynamics(::DecayFactory) = (X=decay_tendency,)

nothing #hide

# ## Construct the model

# We now pass a `DecayFactory()` instance to Agate's generic constructor. The factory type
# selects all the methods defined above, so `construct_factory` can assemble the parameters
# and tracer tendency without containing any Decay-specific code itself.
#
# Agate models normally use `PAR` as an auxiliary field. This model does not need any
# auxiliary fields, so `auxiliary_fields=()` disables that default.

bgc = construct_factory(DecayFactory(); auxiliary_fields=())

# The constructed model requires only the tracer `X`, and the parameter directory contains
# the default decay rate declared above.

println(required_biogeochemical_tracers(bgc))
println("decay rate = ", bgc.parameters.decay_rate * day, " day⁻¹")

# `Val(:X)` uses the same dispatch idea as `Val(:Decay)` above. Here it identifies which
# tracer tendency we want Agate to evaluate. With `X = 2`, the tendency is `-kX`, so for the
# default `k = 0.1 day⁻¹` the result is `-0.2 day⁻¹`.

println("tendency at X = 2: ", bgc(Val(:X), 0, 0, 0, 0, 2.0) * day, " day⁻¹")

nothing #hide