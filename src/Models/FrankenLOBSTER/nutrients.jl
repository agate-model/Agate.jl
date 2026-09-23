"""Source-specific share of one bounded NO3+NH4, Fe-co-limited growth capacity."""

import ...Processes:
    AbstractFactor, AbstractFormulation, FactorComponent, ParameterSlot,
    authored_parameter_bindings, factor_inputs, factor_value, parameter_slots,
    _canonical_bindings
using ...Library.Nutrients: liebig_minimum, monod_limitation

struct ModifiedMonodNitrogenSource{Source} <: AbstractFormulation end

function ModifiedMonodNitrogenSource(source::Symbol)
    source in (:NO₃, :NH₄) || throw(ArgumentError("nitrogen source must be :NO₃ or :NH₄"))
    return ModifiedMonodNitrogenSource{source}()
end

struct NitrogenIronSourceResponse{F<:ModifiedMonodNitrogenSource} <: AbstractFactor
    formulation::F
    nitrate::Symbol
    ammonium::Symbol
    iron::Symbol
    bindings::NamedTuple
end

function NitrogenIronSourceResponse(
    source::Symbol;
    nitrate::Symbol=:NO₃,
    ammonium::Symbol=:NH₄,
    iron::Symbol=:Fe,
    bindings::NamedTuple=NamedTuple(),
)
    return NitrogenIronSourceResponse(
        ModifiedMonodNitrogenSource(source), nitrate, ammonium, iron,
        _canonical_bindings(bindings),
    )
end

authored_parameter_bindings(factor::NitrogenIronSourceResponse) = factor.bindings
factor_inputs(factor::NitrogenIronSourceResponse) = (
    FactorComponent(factor.nitrate), FactorComponent(factor.ammonium), FactorComponent(factor.iron),
)

parameter_slots(::ModifiedMonodNitrogenSource) = (
    ParameterSlot(:nitrate_half_saturation, (:plankton,); domain=:nonnegative),
    ParameterSlot(:ammonium_half_saturation, (:plankton,); domain=:nonnegative),
    ParameterSlot(:iron_half_saturation, (:plankton,); domain=:nonnegative),
    ParameterSlot(:ammonium_inhibition; domain=:nonnegative),
)

@inline function factor_value(
    ::ModifiedMonodNitrogenSource{Source}, nitrate, ammonium, iron,
    nitrate_half_saturation, ammonium_half_saturation, iron_half_saturation,
    ammonium_inhibition,
) where Source
    nitrate_response = monod_limitation(nitrate, nitrate_half_saturation) *
                       exp(-ammonium_inhibition * ammonium)
    ammonium_response = max(
        zero(ammonium), monod_limitation(ammonium, ammonium_half_saturation)
    )
    response_sum = nitrate_response + ammonium_response
    response_sum > zero(response_sum) || return zero(response_sum)

    nitrogen_limitation = min(one(response_sum), max(zero(response_sum), response_sum))
    limitation = liebig_minimum(
        nitrogen_limitation, monod_limitation(iron, iron_half_saturation)
    )
    source_response = Source === :NO₃ ? nitrate_response : ammonium_response
    return limitation * source_response / response_sum
end
