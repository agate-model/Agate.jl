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
participant role can still provide ecological applicability without becoming a storage axis.
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
        ArgumentError("parameter slot axes currently support at most two dimensions"),
    )
    all(axis -> axis isa Symbol, axes) || throw(
        ArgumentError("parameter slot axes must contain only Symbols"),
    )
    return ParameterSlot(name, axes, qualify)
end

parameter_slots(::AbstractFormulation) = ()
parameter_slots(::Smith) = (
    ParameterSlot(:maximum_rate, (:population,)),
    ParameterSlot(:alpha, (:population,)),
)
parameter_slots(::Geider) = (
    ParameterSlot(:maximum_rate, (:population,)),
    ParameterSlot(:alpha, (:population,)),
    ParameterSlot(:chlorophyll_to_carbon_ratio, (:population,)),
)
parameter_slots(::Monod) = (ParameterSlot(:K, (:population,); qualify=:resource),)
parameter_slots(::NormalizedDroop) = (
    ParameterSlot(:minimum_quota, (:population,)),
    ParameterSlot(:maximum_quota, (:population,)),
)
parameter_slots(::QuotaRegulatedMonod) = (
    ParameterSlot(:maximum_rate, (:population,)),
    ParameterSlot(:K, (:population,)),
    ParameterSlot(:minimum_quota, (:population,)),
    ParameterSlot(:maximum_quota, (:population,)),
    ParameterSlot(:hill, (:population,)),
)
parameter_slots(::Liebig) = ()
parameter_slots(::FrankTNorm) = (ParameterSlot(:sharpness),)
parameter_slots(::Q10) = (
    ParameterSlot(:q10),
    ParameterSlot(:reference_temperature),
)
parameter_slots(::IdealizedGrazing) = (
    ParameterSlot(:maximum_rate, (:consumer,)),
    ParameterSlot(:half_saturation, (:consumer,)),
    ParameterSlot(:palatability, (:consumer, :resource)),
    ParameterSlot(:assimilation, (:consumer, :resource)),
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
parameter_slots(::LinearMortality) = (ParameterSlot(:rate, (:population,); qualify=:population),)
parameter_slots(::QuadraticMortality) = (ParameterSlot(:rate, (:population,); qualify=:population),)
parameter_slots(::LinearRemineralization) = (
    ParameterSlot(:rate; qualify=:source),
)
parameter_slots(products::Products) = length(products.targets) == 1 ? () :
    (ParameterSlot(:fraction; qualify=:product),)
parameter_slots(::FixedStoichiometry) = (ParameterSlot(:ratio; qualify=:currency),)

"""Resolved setup-time mapping from one local parameter slot to model storage.

The scientific-tree location and qualifier identify the local slot during compilation.
`axes` is the authoritative semantic axis signature supplied by the formulation schema.
"""
struct ParameterBinding{P,Q,A}
    process::Symbol
    path::P
    slot::Symbol
    qualifier::Q
    axes::A
    parameter::Symbol
end

_parameter_slot_source(node::Union{AbstractFormulation,AbstractStoichiometry,Products}) = node
_parameter_slot_source(node) = formulation(node)

