"""Source-specific share of one bounded NO3+NH4 growth capacity."""

import ...Processes:
    AbstractFactor, AbstractFormulation, FactorComponent, ParameterSlot,
    authored_parameter_bindings, factor_inputs, factor_value, parameter_slots,
    _canonical_bindings
using ...Library.Nutrients: monod_limitation

struct ModifiedMonodNitrogenSource{Source} <: AbstractFormulation end

function ModifiedMonodNitrogenSource(source::Symbol)
    source in (:NO₃, :NH₄) || throw(ArgumentError("nitrogen source must be :NO₃ or :NH₄"))
    return ModifiedMonodNitrogenSource{source}()
end

struct NitrogenSourceResponse{F<:ModifiedMonodNitrogenSource} <: AbstractFactor
    formulation::F
    nitrate::Symbol
    ammonium::Symbol
    bindings::NamedTuple
end

function NitrogenSourceResponse(
    source::Symbol;
    nitrate::Symbol=:NO₃,
    ammonium::Symbol=:NH₄,
    bindings::NamedTuple=NamedTuple(),
)
    return NitrogenSourceResponse(
        ModifiedMonodNitrogenSource(source), nitrate, ammonium, _canonical_bindings(bindings)
    )
end

authored_parameter_bindings(factor::NitrogenSourceResponse) = factor.bindings
factor_inputs(factor::NitrogenSourceResponse) = (
    FactorComponent(factor.nitrate), FactorComponent(factor.ammonium),
)

parameter_slots(::ModifiedMonodNitrogenSource) = (
    ParameterSlot(:nitrate_half_saturation, (:plankton,); domain=:nonnegative),
    ParameterSlot(:ammonium_half_saturation, (:plankton,); domain=:nonnegative),
    ParameterSlot(:ammonium_inhibition; domain=:nonnegative),
)

@inline function factor_value(
    ::ModifiedMonodNitrogenSource{Source}, nitrate, ammonium,
    nitrate_half_saturation, ammonium_half_saturation, ammonium_inhibition,
) where Source
    nitrate_response = monod_limitation(nitrate, nitrate_half_saturation) *
                       exp(-ammonium_inhibition * ammonium)
    ammonium_response = max(
        zero(ammonium), monod_limitation(ammonium, ammonium_half_saturation)
    )
    response_sum = nitrate_response + ammonium_response
    response_sum > zero(response_sum) || return zero(response_sum)

    nitrogen_limitation = min(one(response_sum), max(zero(response_sum), response_sum))
    source_response = Source === :NO₃ ? nitrate_response : ammonium_response
    return nitrogen_limitation * source_response / response_sum
end
