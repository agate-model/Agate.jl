"""Lean static runtime IR used by compiled process tendencies."""

"""Static operand that reads one pre-indexed tracer or auxiliary input."""
struct InputOp{Index} end

"""Static operand that reads one pre-indexed runtime parameter value."""
struct ParameterOp{Name,Indices} end

"""Static operand whose value is one minus the sum of child operands."""
struct ComplementOp{Operands}
    operands::Operands
end

"""Static operand that evaluates a tuple of child operands."""
struct TupleOp{Operands}
    operands::Operands
end

@inline operand_value(::InputOp{Index}, bgc, args) where {Index} = @inbounds args[Index]
@inline operand_value(::ParameterOp{Name,()}, bgc, args) where {Name} = getproperty(bgc.parameters, Name)
@inline operand_value(::ParameterOp{Name,Indices}, bgc, args) where {Name,Indices} =
    @inbounds getproperty(bgc.parameters, Name)[Indices...]

@inline function _operand_sum(operands::Tuple{T1}, bgc, args) where {T1}
    return operand_value(first(operands), bgc, args)
end

@inline function _operand_sum(operands::Tuple{T1,T2,Vararg{Any,N}}, bgc, args) where {T1,T2,N}
    return operand_value(first(operands), bgc, args) + _operand_sum(Base.tail(operands), bgc, args)
end

@inline function operand_value(op::ComplementOp, bgc, args)
    total = _operand_sum(op.operands, bgc, args)
    return one(total) - total
end
@inline operand_value(op::TupleOp, bgc, args) = operand_values(op.operands, bgc, args)

@inline operand_values(::Tuple{}, bgc, args) = ()
@inline function operand_values(operands::Tuple, bgc, args)
    return (operand_value(first(operands), bgc, args), operand_values(Base.tail(operands), bgc, args)...)
end

"""One multiplicative factor attached to a canonical process rate."""
struct FactorOp{Formulation,Operands}
    formulation::Formulation
    operands::Operands
end

@inline function factor_op_value(factor::FactorOp, bgc, args)
    values = operand_values(factor.operands, bgc, args)
    return factor_value(factor.formulation, values...)
end

@inline operand_value(factor::FactorOp, bgc, args) =
    factor_op_value(factor, bgc, args)

@inline apply_rate_factors(::Tuple{}, bgc, args, rate) = rate
@inline function apply_rate_factors(factors::Tuple, bgc, args, rate)
    value = factor_op_value(first(factors), bgc, args)
    return apply_rate_factors(Base.tail(factors), bgc, args, rate * value)
end

"""Canonical runtime process rate: formulation plus static operand tuple."""
struct RateOp{Formulation,Operands,Factors}
    formulation::Formulation
    operands::Operands
    factors::Factors
end

RateOp(formulation, operands::Tuple; factors=()) =
    RateOp{typeof(formulation),typeof(operands),typeof(factors)}(formulation, operands, factors)

@inline function (rate::RateOp)(bgc, args)
    values = operand_values(rate.operands, bgc, args)
    value = process_rate(rate.formulation, values...)
    return apply_rate_factors(rate.factors, bgc, args, value)
end

"""Static multiplicative flux weight with sign encoded in the type."""
struct Weight{Sign,Operands}
    operands::Operands
end

Weight{Sign}(operands::Tuple=()) where {Sign} = Weight{Sign,typeof(operands)}(operands)

@inline function weight_product(operands::Tuple{T1}, bgc, args) where {T1}
    return operand_value(first(operands), bgc, args)
end

@inline function weight_product(operands::Tuple{T1,T2,Vararg{Any,N}}, bgc, args) where {T1,T2,N}
    return operand_value(first(operands), bgc, args) * weight_product(Base.tail(operands), bgc, args)
end

@inline apply_weight(::Weight{1,Tuple{}}, bgc, args, rate) = rate
@inline apply_weight(::Weight{-1,Tuple{}}, bgc, args, rate) = -rate
@inline apply_weight(weight::Weight{1}, bgc, args, rate) = weight_product(weight.operands, bgc, args) * rate
@inline apply_weight(weight::Weight{-1}, bgc, args, rate) = -weight_product(weight.operands, bgc, args) * rate

"""Lean runtime flux term containing only its scientific rate and weight."""
struct FluxTerm{Rate,FluxWeight}
    rate::Rate
    weight::FluxWeight
end

@inline function (term::FluxTerm)(bgc, args)
    rate = term.rate(bgc, args)
    return apply_weight(term.weight, bgc, args, rate)
end

struct StaticFluxEquation{Terms}
    terms::Terms
end

struct StaticZeroEquation end

@inline (::StaticZeroEquation)(bgc, x, y, z, t, args...) = zero(first(args))

@inline function _sum_fluxes(terms::Tuple{T1}, bgc, args) where {T1}
    return first(terms)(bgc, args)
end

@inline function _sum_fluxes(terms::Tuple{T1,T2,Vararg{Any,N}}, bgc, args) where {T1,T2,N}
    return first(terms)(bgc, args) + _sum_fluxes(Base.tail(terms), bgc, args)
end

@inline function (equation::StaticFluxEquation)(bgc, x, y, z, t, args...)
    return _sum_fluxes(equation.terms, bgc, args)
end

