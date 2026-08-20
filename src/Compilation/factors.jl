function _temperature_factor(named::NamedProcess)
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
    return name => factor
end

function _rate_factors(definition::NormalizedModelDefinition, named::NamedProcess)
    match = _temperature_factor(named)
    isnothing(match) && return ()
    name, factor = match
    slots = parameter_slot_bindings(
        definition, named, (:factors, name), factor.formulation
    )
    operands = (
        TracerOp{factor.driver}(),
        parameter_operand(slots.q10),
        parameter_operand(slots.reference_temperature),
    )
    return (FactorElement(factor.formulation, operands),)
end
