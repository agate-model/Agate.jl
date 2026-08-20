"""Resolved setup-time binding for one Q10 temperature factor."""
struct TemperatureFactorBinding{F}
    formulation::F
    driver::Symbol
    q10::Symbol
    reference_temperature::Symbol
end

function _temperature_factor_binding(
    definition::NormalizedModelDefinition, named::NamedProcess
)
    matches = ()
    for pair in pairs(factors(named))
        last(pair) isa Temperature && (matches = (matches..., pair))
    end
    isempty(matches) && return nothing
    length(matches) == 1 || throw(
        ArgumentError("process :$(process_id(named)) must declare at most one temperature factor"),
    )
    name, factor = only(matches)
    factor isa Temperature{Q10} || throw(
        ArgumentError("unsupported temperature factor $(typeof(factor))"),
    )
    path = (:factors, name)
    q10 = parameter_name(definition, _parameter_requirement(named, path, :q10))
    reference = parameter_name(
        definition, _parameter_requirement(named, path, :reference_temperature)
    )
    return TemperatureFactorBinding(factor.formulation, factor.driver, q10, reference)
end

_rate_factors(::Nothing) = ()
function _rate_factors(binding::TemperatureFactorBinding)
    operands = (
        TracerOp{binding.driver}(),
        ScalarParamOp{binding.q10}(),
        ScalarParamOp{binding.reference_temperature}(),
    )
    return (FactorElement(binding.formulation, operands),)
end
