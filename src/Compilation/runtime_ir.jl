"""Lean static runtime IR used by compiled process tendencies."""

"""Static operand that reads one pre-indexed tracer or auxiliary input."""
struct InputOp{I} end

"""Static operand that reads one pre-indexed runtime parameter value."""
struct ParameterOp{S,I} end

"""Static operand whose value is one minus the sum of child operands."""
struct ComplementOp{O}
    operands::O
end

"""Static operand that evaluates a tuple of child operands."""
struct TupleOp{O}
    operands::O
end

@inline operand_value(::InputOp{I}, bgc, args) where {I} = @inbounds args[I]
@inline operand_value(::ParameterOp{S,()}, bgc, args) where {S} = getproperty(bgc.parameters, S)
@inline operand_value(::ParameterOp{S,I}, bgc, args) where {S,I} =
    @inbounds getproperty(bgc.parameters, S)[I...]

@inline function _operand_sum(operands::Tuple{T}, bgc, args) where {T}
    return operand_value(first(operands), bgc, args)
end

@inline function _operand_sum(operands::Tuple{T,S,Vararg{Any,N}}, bgc, args) where {T,S,N}
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
struct FactorOp{F,O}
    formulation::F
    operands::O
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
struct RateOp{F,O,X}
    formulation::F
    operands::O
    factors::X
end

RateOp(formulation, operands::Tuple; factors=()) =
    RateOp{typeof(formulation),typeof(operands),typeof(factors)}(formulation, operands, factors)

@inline function (rate::RateOp)(bgc, args)
    values = operand_values(rate.operands, bgc, args)
    value = process_rate(rate.formulation, values...)
    return apply_rate_factors(rate.factors, bgc, args, value)
end

"""Static multiplicative flux weight with sign encoded in the type."""
struct Weight{Sign,O}
    operands::O
end

Weight{Sign}(operands::Tuple=()) where {Sign} = Weight{Sign,typeof(operands)}(operands)

@inline function weight_product(operands::Tuple{T}, bgc, args) where {T}
    return operand_value(first(operands), bgc, args)
end

@inline function weight_product(operands::Tuple{T,S,Vararg{Any,N}}, bgc, args) where {T,S,N}
    return operand_value(first(operands), bgc, args) * weight_product(Base.tail(operands), bgc, args)
end

@inline apply_weight(::Weight{1,Tuple{}}, bgc, args, rate) = rate
@inline apply_weight(::Weight{-1,Tuple{}}, bgc, args, rate) = -rate
@inline apply_weight(weight::Weight{1}, bgc, args, rate) = weight_product(weight.operands, bgc, args) * rate
@inline apply_weight(weight::Weight{-1}, bgc, args, rate) = -weight_product(weight.operands, bgc, args) * rate

"""Lean runtime flux term containing only its scientific rate and weight."""
struct FluxTerm{R,W}
    rate::R
    weight::W
end

@inline function (term::FluxTerm)(bgc, args)
    rate = term.rate(bgc, args)
    return apply_weight(term.weight, bgc, args, rate)
end

struct StaticFluxEquation{T}
    terms::T
end

struct StaticZeroEquation end

@inline (::StaticZeroEquation)(bgc, x, y, z, t, args...) = zero(first(args))

@inline function _sum_fluxes(terms::Tuple{T}, bgc, args) where {T}
    return first(terms)(bgc, args)
end

@inline function _sum_fluxes(terms::Tuple{T,S,Vararg{Any,N}}, bgc, args) where {T,S,N}
    return first(terms)(bgc, args) + _sum_fluxes(Base.tail(terms), bgc, args)
end

@inline function (equation::StaticFluxEquation)(bgc, x, y, z, t, args...)
    return _sum_fluxes(equation.terms, bgc, args)
end

