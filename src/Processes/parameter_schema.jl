"""Scientific parameter-slot schema for formulations and product mappings."""

"""Zero-or-one setup-time qualifier attached to a local parameter slot."""
struct Qualifier
    axis::Symbol
    value::Symbol
end

"""Formulation-local declaration of one semantic parameter slot.

Dimensionality is structural: zero, one, or two declared axes imply scalar, vector, or
matrix values respectively. `qualify` identifies repeated semantic instances without
changing storage dimensionality. `domain` owns the formulation's numeric validity contract.
For a scalar slot, a qualifier that is also a process participant role provides ecological
applicability without becoming a storage axis.
"""
struct ParameterSlot{Axes<:Tuple}
    name::Symbol
    axes::Axes
    qualify::Union{Nothing,Symbol}
    domain::Symbol
end

function ParameterSlot(
    name::Symbol,
    axes::Tuple=();
    qualify::Union{Nothing,Symbol}=nothing,
    domain::Symbol=:finite,
)
    length(axes) <= 2 || throw(
        ArgumentError("parameter slot axes support at most two dimensions"),
    )
    all(axis -> axis isa Symbol, axes) || throw(
        ArgumentError("parameter slot axes must contain only Symbols"),
    )
    allunique(axes) || throw(ArgumentError("parameter slot axes must be unique"))
    domain in (:finite, :nonnegative, :positive, :unit_interval) || throw(
        ArgumentError("parameter slot domain must be :finite, :nonnegative, :positive, or :unit_interval"),
    )
    return ParameterSlot(name, axes, qualify, domain)
end

"""Return semantic parameter slots declared by a slot-owning scientific node."""
parameter_slots(::AbstractFormulation) = ()
parameter_slots(::FactorizedGrowth) = (
    ParameterSlot(:maximum_rate, (:plankton,); domain=:nonnegative),
)
parameter_slots(::Smith) = (ParameterSlot(:alpha, (:plankton,); domain=:nonnegative),)
parameter_slots(::Geider) = (
    ParameterSlot(:alpha, (:plankton,); domain=:nonnegative),
    ParameterSlot(:chlorophyll_to_carbon_ratio, (:plankton,); domain=:nonnegative),
)
parameter_slots(::Monod) = (
    ParameterSlot(:half_saturation, (:plankton,); domain=:nonnegative),
)
parameter_slots(::NormalizedDroop) = (
    ParameterSlot(:minimum_quota, (:plankton,); domain=:positive),
    ParameterSlot(:maximum_quota, (:plankton,); domain=:positive),
)
parameter_slots(::QuotaRegulatedMonod) = (
    ParameterSlot(:maximum_rate, (:plankton,); domain=:nonnegative),
    ParameterSlot(:half_saturation, (:plankton,); domain=:nonnegative),
    ParameterSlot(:minimum_quota, (:plankton,); domain=:positive),
    ParameterSlot(:maximum_quota, (:plankton,); domain=:positive),
    ParameterSlot(:hill, (:plankton,); domain=:positive),
)
parameter_slots(::Liebig) = ()
parameter_slots(::FrankTNorm) = (ParameterSlot(:sharpness; domain=:positive),)
parameter_slots(::Q10) = (
    ParameterSlot(:q10; domain=:positive),
    ParameterSlot(:reference_temperature),
)
parameter_slots(::PreferentialGrazing) = (
    ParameterSlot(:maximum_rate, (:consumer,); domain=:nonnegative),
    ParameterSlot(:half_saturation, (:consumer,); domain=:nonnegative),
    ParameterSlot(:palatability, (:consumer, :resource); domain=:nonnegative),
    ParameterSlot(:assimilation, (:consumer, :resource); domain=:unit_interval),
)
parameter_slots(::HeterotrophicConsumption) = (
    ParameterSlot(:maximum_rate, (:consumer,); domain=:nonnegative),
    ParameterSlot(:half_saturation, (:resource,); domain=:nonnegative),
    ParameterSlot(:assimilation, (:consumer, :resource); domain=:unit_interval),
)
parameter_slots(::LinearMortality) = (
    ParameterSlot(:rate, (:plankton,); qualify=:plankton, domain=:nonnegative),
)
parameter_slots(::QuadraticMortality) = (
    ParameterSlot(:rate, (:plankton,); qualify=:plankton, domain=:nonnegative),
)
parameter_slots(::LinearRemineralization) = (
    ParameterSlot(:rate; qualify=:source, domain=:nonnegative),
)
parameter_slots(products::Products) = length(products.destinations) == 1 ? () :
    (ParameterSlot(:fraction; qualify=:product, domain=:unit_interval),)
parameter_slots(::FixedStoichiometry) = (
    ParameterSlot(:ratio; qualify=:element, domain=:nonnegative),
)

"""Resolved setup-time mapping from one local parameter slot to model storage.

`axes` is the authoritative semantic axis signature supplied by the formulation schema, while
`axis_components` records the already-resolved participant components that determine realized
ecological applicability.
"""
struct ParameterBinding{Axes,AxisComponents}
    process::Symbol
    path::Tuple
    slot::Symbol
    axes::Axes
    axis_components::AxisComponents
    parameter::Symbol
    domain::Symbol
end

_parameter_slot_source(node::Union{AbstractFormulation,AbstractStoichiometry,Products}) = node
_parameter_slot_source(node) = formulation(node)

