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

struct UnitFactorTerm end

struct TemperatureFactorTerm{F,D,Q,R}
    formulation::F
end

_lower_factor(::Nothing) = UnitFactorTerm()
function _lower_factor(binding::TemperatureFactorBinding{F}) where {F}
    return TemperatureFactorTerm{
        F,binding.driver,binding.q10,binding.reference_temperature
    }(binding.formulation)
end

@inline _apply_factor(::UnitFactorTerm, bgc, args, rate) = rate

@inline function _apply_factor(
    term::TemperatureFactorTerm{F,D,Q,R}, bgc, args, rate
) where {F,D,Q,R}
    temperature = getproperty(bgc.tracers, D)(args)
    q10 = getproperty(bgc.parameters, Q)
    reference = getproperty(bgc.parameters, R)
    return rate * factor_value(term.formulation, temperature, q10, reference)
end
