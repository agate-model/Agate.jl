using ...Processes:
    AbstractFactor,
    AbstractFormulation,
    FactorDriver,
    FactorComponent,
    ParameterSlot

import ...Processes:
    authored_parameter_bindings,
    factor_inputs,
    factor_value,
    formulation,
    parameter_slots

"""LOBSTER-style saturating PAR response, ``PAR / (K_PAR + PAR)``."""
struct SaturatingLightFormulation <: AbstractFormulation end

struct SaturatingLight <: AbstractFactor
    driver::Symbol
    bindings::NamedTuple
end

SaturatingLight(; driver::Symbol=:PAR, bindings=(half_saturation=:light_half_saturation,)) =
    SaturatingLight(driver, bindings)

formulation(::SaturatingLight) = SaturatingLightFormulation()
authored_parameter_bindings(factor::SaturatingLight) = factor.bindings
parameter_slots(::SaturatingLightFormulation) =
    (ParameterSlot(:half_saturation, (:plankton,); domain=:nonnegative),)
factor_inputs(factor::SaturatingLight) = (FactorDriver(factor.driver),)

@inline function factor_value(::SaturatingLightFormulation, PAR, half_saturation)
    PAR == zero(PAR) && half_saturation == zero(half_saturation) && return zero(PAR)
    return PAR / (half_saturation + PAR)
end

"""Suppression of nitrate uptake by ammonium, ``exp(-psi * NH4)``."""
struct AmmoniumInhibitionFormulation <: AbstractFormulation end

struct AmmoniumInhibition <: AbstractFactor
    resource::Symbol
    bindings::NamedTuple
end

AmmoniumInhibition(;
    resource::Symbol=:NH₄,
    bindings=(coefficient=:nitrate_ammonia_inhibition,),
) = AmmoniumInhibition(resource, bindings)

formulation(::AmmoniumInhibition) = AmmoniumInhibitionFormulation()
authored_parameter_bindings(factor::AmmoniumInhibition) = factor.bindings
parameter_slots(::AmmoniumInhibitionFormulation) =
    (ParameterSlot(:coefficient; domain=:nonnegative),)
factor_inputs(factor::AmmoniumInhibition) = (FactorComponent(factor.resource),)

@inline factor_value(::AmmoniumInhibitionFormulation, ammonium, coefficient) =
    exp(-coefficient * ammonium)
