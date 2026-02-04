"""Citation2026 style DARWIN configuration.

This is an *example* variant showing how to register a paper-specific recipe
without forking the whole DARWIN module.

The intent is to keep it as a thin construction-time wrapper:
- pick defaults (community sizes/diameters)
- use the family defaults for dynamics
- leave all biology implementation in `Models/DARWIN/*`
"""

using ...Interface:
    default_plankton_dynamics, default_biogeochem_dynamics, default_community
using ...Constructor: build_plankton_community
using Agate.Models: ModelId, VariantSpec, register_variant

"""Return the `DARWIN/citation2026/A` variant specification."""
function citation2026_A_spec(;
    n_phyto::Int=20,
    n_zoo::Int=10,
    phyto_diameters=(2, 10, :log_splitting),
    zoo_diameters=(20, 100, :linear_splitting),
    parameters::NamedTuple=(;),
    interaction_overrides::Union{Nothing,NamedTuple}=nothing,
)
    factory = DarwinFactory()

    base = default_community(factory)
    community = build_plankton_community(
        base;
        n=(Z=n_zoo, P=n_phyto),
        diameters=(Z=zoo_diameters, P=phyto_diameters),
    )

    roles = (consumers=(:Z,), prey=(:P,))

    return VariantSpec(
        ModelId(:DARWIN, :citation2026, :A),
        factory,
        default_plankton_dynamics(factory),
        default_biogeochem_dynamics(factory),
        community,
        roles,
        parameters,
        interaction_overrides,
    )
end

register_variant(ModelId(:DARWIN, :citation2026, :A), citation2026_A_spec)
