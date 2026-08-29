"""Scientific parameter-slot schema for formulations and product mappings."""

"""Zero-or-one setup-time qualifier attached to a local parameter slot."""
struct Qualifier
    axis::Symbol
    value::Symbol
end

"""Formulation-local declaration of one semantic parameter slot.

Dimensionality is structural: zero, one, or two declared axes imply scalar, vector, or
matrix values respectively. `qualify` identifies repeated semantic instances without
changing storage dimensionality. For a scalar slot, a qualifier that is also a process
participant role provides ecological applicability without becoming a storage axis.
"""
struct ParameterSlot{A<:Tuple}
    name::Symbol
    axes::A
    qualify::Union{Nothing,Symbol}
end

function ParameterSlot(
    name::Symbol,
    axes::Tuple=();
    qualify::Union{Nothing,Symbol}=nothing,
)
    length(axes) <= 2 || throw(
        ArgumentError("parameter slot axes support at most two dimensions"),
    )
    all(axis -> axis isa Symbol, axes) || throw(
        ArgumentError("parameter slot axes must contain only Symbols"),
    )
    return ParameterSlot(name, axes, qualify)
end

"""Return semantic parameter slots declared by a slot-owning scientific node."""
parameter_slots(::AbstractFormulation) = ()
parameter_slots(::FactorizedGrowth) = (ParameterSlot(:maximum_rate, (:plankton,)),)
parameter_slots(::Smith) = (ParameterSlot(:alpha, (:plankton,)),)
parameter_slots(::Geider) = (
    ParameterSlot(:alpha, (:plankton,)),
    ParameterSlot(:chlorophyll_to_carbon_ratio, (:plankton,)),
)
parameter_slots(::Monod) = (ParameterSlot(:half_saturation, (:plankton,)),)
parameter_slots(::NormalizedDroop) = (
    ParameterSlot(:minimum_quota, (:plankton,)),
    ParameterSlot(:maximum_quota, (:plankton,)),
)
parameter_slots(::QuotaRegulatedMonod) = (
    ParameterSlot(:maximum_rate, (:plankton,)),
    ParameterSlot(:half_saturation, (:plankton,)),
    ParameterSlot(:minimum_quota, (:plankton,)),
    ParameterSlot(:maximum_quota, (:plankton,)),
    ParameterSlot(:hill, (:plankton,)),
)
parameter_slots(::Liebig) = ()
parameter_slots(::FrankTNorm) = (ParameterSlot(:sharpness),)
parameter_slots(::Q10) = (
    ParameterSlot(:q10),
    ParameterSlot(:reference_temperature),
)
parameter_slots(::PreferentialGrazing) = (
    ParameterSlot(:maximum_rate, (:consumer,)),
    ParameterSlot(:half_saturation, (:consumer,)),
    ParameterSlot(:palatability, (:consumer, :resource)),
    ParameterSlot(:assimilation, (:consumer, :resource)),
)
parameter_slots(::HeterotrophicConsumption) = (
    ParameterSlot(:maximum_rate, (:consumer,)),
    ParameterSlot(:half_saturation, (:resource,)),
    ParameterSlot(:assimilation, (:consumer, :resource)),
)
parameter_slots(::LinearMortality) = (ParameterSlot(:rate, (:plankton,); qualify=:plankton),)
parameter_slots(::QuadraticMortality) = (ParameterSlot(:rate, (:plankton,); qualify=:plankton),)
parameter_slots(::LinearRemineralization) = (
    ParameterSlot(:rate; qualify=:source),
)
parameter_slots(products::Products) = length(products.targets) == 1 ? () :
    (ParameterSlot(:fraction; qualify=:product),)
parameter_slots(::FixedStoichiometry) = (ParameterSlot(:ratio; qualify=:element),)

"""Resolved setup-time mapping from one local parameter slot to model storage.

`axes` is the authoritative semantic axis signature supplied by the formulation schema, while
`axis_components` records the already-resolved participant components that determine realized
ecological applicability.
"""
struct ParameterBinding{A,C}
    axes::A
    axis_components::C
    parameter::Symbol
end

_parameter_slot_source(node::Union{AbstractFormulation,AbstractStoichiometry,Products}) = node
_parameter_slot_source(node) = formulation(node)

