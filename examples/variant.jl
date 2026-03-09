# # [DARWIN variant example] (@id DARWIN_variant_example)

# In this example we define a variation of the [Agate.jl-DARWIN](@ref DARWIN) model.
# Variants are a powerful way to implement paper-specific constructors without forking the whole model.
# In this example, we implement a variant of the DARWIN model with a different set of interaction roles.

# ## Loading dependencies
# The example uses the Agate.jl modules Factories (defaults), Configuration (community building), and Models (variant registration).

using Agate.Factories:
    default_plankton_dynamics, default_biogeochem_dynamics, default_community
using Agate.Configuration: build_plankton_community
using Agate.Models: ModelId, VariantSpec, register_variant, variant, construct
using Agate.Models.DARWIN: DarwinFactory

# ## ModelId definition
# First we define a `ModelId` for our variant. The `ModelId` is a stable identifier that uniquely identifies the variant. It consists of three parts:
# - `family`: the model family (e.g. DARWIN, NiPiZD, etc.)
# - `citation`: a citation key that identifies the paper or source of the variant (e.g. citation2026)
# - `tag`: a stable label for the variant within that citation (e.g. A, B, submission, accepted). 

model_id = ModelId(:DARWIN, :citation2026, :B)
nothing #hide

# ## Variant specification
# Next we define a function that returns the variant specification. The specification includes:
# - the model ID (family, paper, variant)
# - the factory to use for construction
# - the dynamics to use for plankton and biogeochemistry
# - the plankton community to use
# - the interaction roles (e.g. which groups are consumers and which are prey)
# - any auxiliary fields to include in the model (e.g. PAR)
# - any parameters to override from the defaults
# In this case, we will only update the interaction_roles, leaving the rest as defaults.

function citation2026_B_spec(;
    n_phyto::Int=2,
    n_zoo::Int=2,
    phyto_diameters=(2, 10, :log_splitting),
    zoo_diameters=(20, 100, :linear_splitting),
    parameters::NamedTuple=(;),
    interaction_overrides::Union{Nothing,NamedTuple}=nothing,
    auxiliary_fields::Tuple=(:PAR,),
)
    factory = DarwinFactory()

    base = default_community(factory)
    community = build_plankton_community(
        base; n=(Z=n_zoo, P=n_phyto), diameters=(Z=zoo_diameters, P=phyto_diameters)
    )

    interaction_roles = (consumers=(:Z,), prey=(:P, :Z)) # allow Z to eat Z 

    return VariantSpec(
        model_id,
        factory,
        default_plankton_dynamics(factory),
        default_biogeochem_dynamics(factory),
        community,
        interaction_roles,
        auxiliary_fields,
        parameters,
        interaction_overrides,
    )
end
nothing #hide

# ## Registering the variant
#
# Registering populates an in-memory registry. You can then look up the variant by `ModelId`.
#
# Important:
# - Registration must run **in the current Julia session** before calling `variant(model_id)`.
# - If this lives in `src/Models/DARWIN/Variants/*.jl` and is included by `Variants/variants.jl`,
#   then `using Agate.Models.DARWIN` will register it automatically at module load.

id = register_variant(model_id, citation2026_B_spec)
nothing #hide

# ## Constructing the variant
#
# 1) `variant(id; kwargs...)` returns the `VariantSpec` (kwargs forwarded to the spec builder)
# 2) `construct(spec; ...)` constructs the concrete biogeochemistry model

spec = variant(id; n_phyto=3, n_zoo=2)
bgc = construct(spec)
nothing #hide

# ## Calling existing variants by ModelId
#
# You can always look up by ModelId — as long as the code that called `register_variant(...)`
# has run in this session.
#
# For built-in variants shipped with Agate, you typically just load the model module:

using Agate.Models.DARWIN  # ensures DARWIN’s built-in variants are registered

built_in_id = ModelId(:DARWIN, :citation2026, :A)  # example
built_in_spec = variant(built_in_id)              # looks up and calls the registered builder
built_in_bgc = construct(built_in_spec)
nothing #hide

# ## Calling your own variants
#
# If you define a variant in your own project (not inside Agate’s src tree),
# you must `include` (or `using MyVariantsPackage`) before calling `variant(model_id)`.
