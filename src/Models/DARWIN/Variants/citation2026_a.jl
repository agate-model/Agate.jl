# """Citation2026 style DARWIN configuration.
# 
# This is an *example* variant showing how to register a paper-specific recipe
# without forking the whole DARWIN module.
# 
# The intent is to keep it as a thin construction-time wrapper:
# - pick defaults (community sizes/diameters)
# - use the family defaults for dynamics
# - leave all biology implementation in `Models/DARWIN/*`
# """

using ...Factories:
    default_plankton_dynamics, default_biogeochem_dynamics, default_community
using ...Configuration: build_plankton_community
using Agate.Models: ModelId, VariantSpec, register_variant

"""Return the `DARWIN/citation2026/A` variant specification."""
function citation2026_A_spec(;
    phyto_diameters=(n=20, min_esd=2, max_esd=10, splitting=:log_splitting),
    zoo_diameters=(n=10, min_esd=20, max_esd=100, splitting=:linear_splitting),
    parameters::NamedTuple=(;),
    interaction_overrides::Union{Nothing,NamedTuple}=nothing,
    auxiliary_fields::Tuple=(:PAR,),
)
    factory = DarwinFactory()

    base = default_community(factory)
    community = build_plankton_community(
        base; diameters=(Z=zoo_diameters, P=phyto_diameters)
    )

    interaction_roles = (consumers=(:Z,), prey=(:P,))

    return VariantSpec(
        ModelId(:DARWIN, :citation2026, :A),
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

register_variant(ModelId(:DARWIN, :citation2026, :A), citation2026_A_spec)
