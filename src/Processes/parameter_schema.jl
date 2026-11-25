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
    length(axes) <= 2 ||
        throw(ArgumentError("parameter slot axes support at most two dimensions"))
    all(axis -> axis isa Symbol, axes) ||
        throw(ArgumentError("parameter slot axes must contain only Symbols"))
    allunique(axes) || throw(ArgumentError("parameter slot axes must be unique"))
    domain in (:finite, :nonnegative, :positive, :unit_interval) || throw(
        ArgumentError(
            "parameter slot domain must be :finite, :nonnegative, :positive, or :unit_interval",
        ),
    )
    return ParameterSlot(name, axes, qualify, domain)
end

"""Return semantic parameter slots declared by a slot-owning scientific node."""
parameter_slots(::AbstractFormulation) = ()
function parameter_slots(::FactorizedGrowth)
    return (ParameterSlot(:maximum_rate, (:plankton,); domain=:nonnegative),)
end
function parameter_slots(process::Growth)
    return if isnothing(process.products)
        parameter_slots(FactorizedGrowth())
    else
        (
            parameter_slots(FactorizedGrowth())...,
            ParameterSlot(:product_fraction, (:plankton,); domain=:unit_interval),
        )
    end
end
parameter_slots(::Smith) = (ParameterSlot(:alpha, (:plankton,); domain=:nonnegative),)
function parameter_slots(::Geider)
    return (
        ParameterSlot(:alpha, (:plankton,); domain=:nonnegative),
        ParameterSlot(:chlorophyll_to_carbon_ratio, (:plankton,); domain=:nonnegative),
    )
end
function parameter_slots(::ExponentialSaturation)
    return (ParameterSlot(:light_scale, (:plankton,); domain=:positive),)
end
function parameter_slots(::Monod)
    return (ParameterSlot(:half_saturation, (:plankton,); domain=:nonnegative),)
end
function parameter_slots(::InhibitedMonod)
    return (
        ParameterSlot(:half_saturation, (:plankton,); domain=:nonnegative),
        ParameterSlot(:inhibition; domain=:nonnegative),
    )
end
function parameter_slots(::NormalizedDroop)
    return (
        ParameterSlot(:minimum_quota, (:plankton,); domain=:positive),
        ParameterSlot(:maximum_quota, (:plankton,); domain=:positive),
    )
end
function parameter_slots(::QuotaRegulatedMonod)
    return (
        ParameterSlot(:maximum_rate, (:plankton,); domain=:nonnegative),
        ParameterSlot(:half_saturation, (:plankton,); domain=:nonnegative),
        ParameterSlot(:minimum_quota, (:plankton,); domain=:positive),
        ParameterSlot(:maximum_quota, (:plankton,); domain=:positive),
        ParameterSlot(:hill, (:plankton,); domain=:positive),
    )
end
parameter_slots(::Liebig) = ()
parameter_slots(::FrankTNorm) = (ParameterSlot(:sharpness; domain=:positive),)
function parameter_slots(::Q10{nothing})
    return (ParameterSlot(:q10; domain=:positive), ParameterSlot(:reference_temperature))
end
function parameter_slots(::Q10{Axis}) where {Axis}
    return (
        ParameterSlot(:q10, (Axis,); domain=:positive),
        ParameterSlot(:reference_temperature),
    )
end
function parameter_slots(::PreferentialGrazing)
    return (
        ParameterSlot(:maximum_rate, (:consumer,); domain=:nonnegative),
        ParameterSlot(:half_saturation, (:consumer,); domain=:nonnegative),
        ParameterSlot(:palatability, (:consumer, :resource); domain=:nonnegative),
        ParameterSlot(:assimilation, (:consumer, :resource); domain=:unit_interval),
    )
end
function parameter_slots(::LinearGrazing)
    return (
        ParameterSlot(:rate, (:consumer,); domain=:nonnegative),
        ParameterSlot(:palatability, (:consumer, :resource); domain=:nonnegative),
        ParameterSlot(:assimilation, (:consumer, :resource); domain=:unit_interval),
    )
end
function parameter_slots(::HeterotrophicConsumption)
    return (
        ParameterSlot(:maximum_rate, (:consumer,); domain=:nonnegative),
        ParameterSlot(:half_saturation, (:consumer, :resource); domain=:positive),
        ParameterSlot(:substrate_preference, (:consumer, :resource); domain=:nonnegative),
        ParameterSlot(:assimilation, (:consumer, :resource); domain=:unit_interval),
    )
end
function parameter_slots(::LinearMortality)
    return (ParameterSlot(:rate, (:plankton,); qualify=:plankton, domain=:nonnegative),)
end
function parameter_slots(::QuadraticMortality)
    return (ParameterSlot(:rate, (:plankton,); qualify=:plankton, domain=:nonnegative),)
end
function parameter_slots(::LinearRemineralization)
    return (ParameterSlot(:rate; qualify=:source, domain=:nonnegative),)
end
function parameter_slots(products::Products)
    return if length(products.destinations) == 1
        ()
    else
        (ParameterSlot(:fraction; qualify=:product, domain=:unit_interval),)
    end
end

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

_parameter_slot_source(node::Union{AbstractFormulation,Products,Growth}) = node
_parameter_slot_source(node) = formulation(node)
