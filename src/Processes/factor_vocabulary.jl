"""Abstract supertype for scientific process objects."""
abstract type AbstractProcess end

"""Abstract supertype for concrete scientific formulations."""
abstract type AbstractFormulation end

"""Smith light-limitation formulation."""
struct Smith <: AbstractFormulation end

"""Geider light-response formulation.

This factor regulates Growth rate using chlorophyll-to-carbon information; it does not itself
synthesize a prognostic non-elemental state such as `:chlorophyll`.
"""
struct Geider <: AbstractFormulation end

"""Monod single-resource limitation formulation."""
struct Monod <: AbstractFormulation end

"""Normalized Droop cellular-quota growth-limitation formulation."""
struct NormalizedDroop <: AbstractFormulation end

"""Monod nutrient-uptake formulation regulated by cellular quota capacity."""
struct QuotaRegulatedMonod <: AbstractFormulation end

"""Liebig minimum response-combination formulation."""
struct Liebig <: AbstractFormulation end

"""Differentiable Frank t-norm nutrient-combination formulation."""
struct FrankTNorm <: AbstractFormulation end

"""Q10 temperature-response formulation."""
struct Q10 <: AbstractFormulation end

"""Growth formulation with a base maximum rate and optional multiplicative factors."""
struct FactorizedGrowth <: AbstractFormulation end

"""Abstract supertype for named multiplicative process-rate factors."""
abstract type AbstractFactor end

"""Living-prey grazing formulation with pairwise consumer-resource capacity.

`maximum_rate` is a per-consumer rate applied independently to each declared prey edge, so
each edge receives the full per-consumer rate. `palatability` is a nonnegative interaction weight,
not a probability.
"""
struct PreferentialGrazing <: AbstractFormulation end

"""Heterotrophic resource-consumption formulation with pairwise consumer-resource capacity.

`maximum_rate` is a per-consumer rate applied independently to each declared resource edge, so
each edge receives the full per-consumer rate.
"""
struct HeterotrophicConsumption <: AbstractFormulation end

"""Linear plankton-mortality formulation."""
struct LinearMortality <: AbstractFormulation end

"""Quadratic plankton-mortality formulation."""
struct QuadraticMortality <: AbstractFormulation end

"""Linear source-to-destination remineralization formulation."""
struct LinearRemineralization <: AbstractFormulation end

"""Abstract supertype for process stoichiometry mappings."""
abstract type AbstractStoichiometry end

# Shared authoring canonicalizers

function _canonical_bindings(bindings::NamedTuple)
    names = sort!(collect(keys(bindings)); by=String)
    names_tuple = Tuple(names)
    binding_values = Tuple(begin
        value = getproperty(bindings, name)
        if value isa Symbol
            value
        elseif value isa NamedTuple
            all(entry -> entry isa Symbol, values(value)) || throw(
                ArgumentError(
                    "binding map for :$name must map qualifier values to parameter Symbols"
                ),
            )
            qualifier_names = sort!(collect(keys(value)); by=String)
            qualifier_names_tuple = Tuple(qualifier_names)
            qualifier_values = Tuple(getproperty(value, key) for key in qualifier_names)
            NamedTuple{qualifier_names_tuple}(qualifier_values)
        else
            throw(ArgumentError(
                "binding :$name must be a parameter Symbol or one-level qualifier NamedTuple",
            ))
        end
    end for name in names)
    return NamedTuple{names_tuple}(binding_values)
end

"""Return setup-only authored model-parameter bindings for one slot-owning node."""
function authored_parameter_bindings end

"""Fixed conversion from one reference element to process target elements.

Each bound ratio is the amount of its target element per unit reference element.
"""
struct FixedStoichiometry <: AbstractStoichiometry
    reference_element::Symbol
    bindings::NamedTuple
end

FixedStoichiometry(; reference_element::Symbol, bindings::NamedTuple=NamedTuple()) =
    FixedStoichiometry(reference_element, _canonical_bindings(bindings))

authored_parameter_bindings(stoichiometry::FixedStoichiometry) = stoichiometry.bindings

function _canonical_participants(role::Symbol, values)
    values isa Symbol && (values = (values,))
    values isa Tuple || throw(ArgumentError("participant `$role` must be a Symbol or tuple"))
    isempty(values) && throw(ArgumentError("participant `$role` cannot be empty"))
    all(value -> value isa Symbol, values) || throw(
        ArgumentError("participant `$role` must contain only Symbols"),
    )
    allunique(values) || throw(
        ArgumentError("participant `$role` must not contain duplicates"),
    )
    return values
end

"""Light-dependent multiplicative Growth factor using the Growth rate scale."""
struct Light{Formulation<:Union{Smith,Geider}} <: AbstractFactor
    formulation::Formulation
    driver::Symbol
    bindings::NamedTuple
end

function Light(
    formulation::Union{Smith,Geider}; driver::Symbol, bindings::NamedTuple=NamedTuple()
)
    return Light(formulation, driver, _canonical_bindings(bindings))
end

authored_parameter_bindings(factor::Light) = factor.bindings

"""Single-resource multiplicative nutrient-rate factor.

The factor reads an environmental Pool but does not define process material transfer.
"""
struct NutrientResponse{Formulation<:Monod} <: AbstractFactor
    formulation::Formulation
    resource::Symbol
    bindings::NamedTuple
end

function NutrientResponse(
    formulation::Monod; resource::Symbol, bindings::NamedTuple=NamedTuple()
)
    return NutrientResponse(formulation, resource, _canonical_bindings(bindings))
end

authored_parameter_bindings(factor::NutrientResponse) = factor.bindings

"""Cellular-quota response used by quota-limited growth.

`variable_state` identifies the internal inventory whose quota varies relative to the
Growth plankton's intrinsic reference state.
"""
struct QuotaResponse{Formulation<:NormalizedDroop} <: AbstractFactor
    formulation::Formulation
    variable_state::Symbol
    bindings::NamedTuple
end

function QuotaResponse(
    formulation::NormalizedDroop;
    variable_state::Symbol,
    bindings::NamedTuple=NamedTuple(),
)
    return QuotaResponse(formulation, variable_state, _canonical_bindings(bindings))
end

authored_parameter_bindings(factor::QuotaResponse) = factor.bindings

"""Temperature-dependent multiplicative process-rate factor."""
struct Temperature{Formulation<:Q10} <: AbstractFactor
    formulation::Formulation
    driver::Symbol
    bindings::NamedTuple
end

function Temperature(
    formulation::Q10; driver::Symbol=:temperature, bindings::NamedTuple=NamedTuple()
)
    return Temperature(formulation, driver, _canonical_bindings(bindings))
end

authored_parameter_bindings(factor::Temperature) = factor.bindings

function _canonical_namedtuple(values::NamedTuple)
    names = sort!(collect(keys(values)); by=String)
    return NamedTuple{Tuple(names)}(Tuple(getproperty(values, name) for name in names))
end

"""Multi-response nutrient factor with formulation-owned response composition.

External `NutrientResponse` subfactors read environmental resource Pools, while internal
`QuotaResponse` subfactors read prognostic cellular states. Each `responses` key is the Element
identity represented by that response (for example, `nitrogen=...`). Both modify process rate
only; material transfer is owned by the process itself, so external and internal responses may
be combined within one `NutrientLimitation` factor.
"""
struct NutrientLimitation{
    Formulation<:Union{Liebig,FrankTNorm},Responses<:NamedTuple
} <: AbstractFactor
    formulation::Formulation
    responses::Responses
    bindings::NamedTuple

    function NutrientLimitation(
        formulation::Formulation, responses::Responses, bindings::NamedTuple
    ) where {Formulation<:Union{Liebig,FrankTNorm},Responses<:NamedTuple}
        isempty(responses) && throw(ArgumentError("nutrient `responses` cannot be empty"))
        all(
            response -> response isa Union{NutrientResponse,QuotaResponse}, values(responses)
        ) || throw(ArgumentError(
            "nutrient `responses` values must be NutrientResponse or QuotaResponse factors",
        ))
        return new{Formulation,Responses}(formulation, responses, bindings)
    end
end

function NutrientLimitation(
    formulation::Union{Liebig,FrankTNorm};
    responses::NamedTuple,
    bindings::NamedTuple=NamedTuple(),
)
    return NutrientLimitation(
        formulation, _canonical_namedtuple(responses), _canonical_bindings(bindings)
    )
end

authored_parameter_bindings(factor::NutrientLimitation) = factor.bindings

abstract type AbstractFactorInput end

"""External driver read required by one scientific factor."""
struct FactorDriver <: AbstractFactorInput
    identity::Symbol
end

"""Scalar model-component read required by one scientific factor."""
struct FactorComponent <: AbstractFactorInput
    component::Symbol
end

"""Setup-only read of one prognostic state from the factor's current logical plankton.

`reference.plankton` must match the logical plankton at the realized `:plankton` axis position.
"""
struct FactorPlanktonState <: AbstractFactorInput
    reference::PlanktonStateRef
end

"""Return the ordered semantic inputs read by a factor before its parameter slots."""
factor_inputs(::AbstractFactor) = ()
factor_inputs(factor::Light) = (FactorDriver(factor.driver),)
factor_inputs(factor::Temperature) = (FactorDriver(factor.driver),)
factor_inputs(factor::NutrientResponse) = (FactorComponent(factor.resource),)
factor_inputs(::QuotaResponse) = ()

"""Return named child factors composed by a factor."""
factor_subfactors(::AbstractFactor) = NamedTuple()
factor_subfactors(factor::NutrientLimitation) = factor.responses

factor_subfactor_path(path::Tuple, ::AbstractFactor, name::Symbol) = (path..., name)
factor_subfactor_path(path::Tuple, ::NutrientLimitation, name::Symbol) = (path..., :responses, name)
