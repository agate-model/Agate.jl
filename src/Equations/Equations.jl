"""Agate.Equations

Construction-time symbolic building blocks for authoring biogeochemical dynamics
without hand-written `Expr` quoting and without any `@register_*` macros.

This module is intentionally non-circular and does not depend on Models or
Constructor. Objects defined here are used *only at construction time*.

Design notes
------------
- Dynamics authors write equations as functions that accept a `PV` placeholder
  namespace (passed in by `construct`) and then use `PV.<name>` references like
  `PV.maximum_growth_rate[i]` or `PV.detritus_remineralization`.

- Those `PV.<name>` bindings are small `ParamVar` placeholders.
- When a `ParamVar` is indexed, it produces an `AExpr` node and records a
  requirement (scalar/vector/matrix) for the parameter key.
- At runtime, parameters are plain numeric scalars and arrays; the symbolic
  placeholders never enter kernels.

Public API:
- `sum_over(items) do sym, idx ... end` — construction-time symbolic sum builder.
- `Equation(::AExpr)` — wrapper returned by dynamics builders.

The old `group_param`/`community_param` distinction is intentionally removed;
missing/`nothing` policy is specified in the parameter registry.
"""

module Equations

using Logging

export Requirements, req, merge_requirements
export AExpr, Equation, expr, requirements
export ParamVar
export sum_over

# -----------------------------------------------------------------------------
# Requirements
# -----------------------------------------------------------------------------

"""Construction-time dependency set for a symbolic expression."""
struct Requirements
    vectors::Vector{Symbol}
    matrices::Vector{Symbol}
    scalars::Vector{Symbol}
end

"""Create a `Requirements` object."""
function req(; vectors=(), matrices=(), scalars=())
    return Requirements(
        Symbol[vectors...],
        Symbol[matrices...],
        Symbol[scalars...],
    )
end

@inline function _append_unique!(dest::Vector{Symbol}, src::Vector{Symbol})
    for s in src
        s in dest || push!(dest, s)
    end
    return dest
end

"""Merge two `Requirements` with stable ordering (first-seen wins)."""
function merge_requirements(r1::Requirements, r2::Requirements)
    return Requirements(
        _append_unique!(copy(r1.vectors), r2.vectors),
        _append_unique!(copy(r1.matrices), r2.matrices),
        _append_unique!(copy(r1.scalars), r2.scalars),
    )
end

# -----------------------------------------------------------------------------
# AExpr (AST node + requirements)
# -----------------------------------------------------------------------------

"""Symbolic expression node plus its construction-time requirements."""
struct AExpr
    node::Any
    req::Requirements
end

@inline AExpr(node) = AExpr(node, req())


"""
    sum_expr(terms)

Return an expression that sums a collection of AST terms.

Terms may be `Expr`, `Symbol`, or literal values (e.g. numbers).
If `terms` is empty, returns `0`.
"""
function sum_expr(terms::AbstractVector)
    isempty(terms) && return 0

    s = terms[1]
    for i in 2:length(terms)
        s = :($s + $(terms[i]))
    end
    return s
end

# -----------------------------------------------------------------------------
# Parameter placeholder
# -----------------------------------------------------------------------------

"""A construction-time placeholder for a parameter named `K`.

A `ParamVar` behaves like a scalar/array reference during equation authoring:
- using it as a scalar records `K` in `requirements(...).scalars`
- indexing with one index records `K` in `requirements(...).vectors`
- indexing with two indices records `K` in `requirements(...).matrices`

At runtime, the generated tracer methods bind `K` to a numeric value/array.
"""
struct ParamVar{K} end

@inline function _to_aexpr(x)
    return x isa AExpr ? x : AExpr(x, req())
end

@inline function _to_aexpr(::ParamVar{K}) where {K}
    return AExpr(K, req(scalars=(K,)))
end

function Base.getindex(::ParamVar{K}, i) where {K}
    return AExpr(Expr(:ref, K, i), req(vectors=(K,)))
end

function Base.getindex(::ParamVar{K}, j, i) where {K}
    return AExpr(Expr(:ref, K, j, i), req(matrices=(K,)))
end

# -----------------------------------------------------------------------------
# Operator overloading
# -----------------------------------------------------------------------------

@inline function _binop(op::Symbol, a, b)
    aa = _to_aexpr(a)
    bb = _to_aexpr(b)
    return AExpr(Expr(:call, op, aa.node, bb.node), merge_requirements(aa.req, bb.req))
end

@inline function _unop(op::Symbol, a)
    aa = _to_aexpr(a)
    return AExpr(Expr(:call, op, aa.node), aa.req)
end

# Minimal arithmetic surface: accept only expression-like inputs.
#
# These methods are deliberately narrow so they do not interfere with ordinary
# numeric arithmetic. When both operands are numbers, Julia dispatches to Base.
const _EquationExprLike = Union{AExpr, ParamVar, Symbol, Number}

Base.:+(a::_EquationExprLike, b::_EquationExprLike) = _binop(:+, a, b)
Base.:-(a::_EquationExprLike, b::_EquationExprLike) = _binop(:-, a, b)
Base.:*(a::_EquationExprLike, b::_EquationExprLike) = _binop(:*, a, b)
Base.:/(a::_EquationExprLike, b::_EquationExprLike) = _binop(:/, a, b)
Base.:^(a::_EquationExprLike, b::_EquationExprLike) = _binop(:^, a, b)
Base.:-(a::_EquationExprLike) = _unop(:-, a)

# -----------------------------------------------------------------------------
# sum_over sum builder (construction-time only)
# -----------------------------------------------------------------------------

"""
    sum_over(items) do sym, idx
        ... -> AExpr
    end

Build a symbolic unrolled sum over `items` (construction-time only).
"""
function sum_over(items, f::Function)
    nodes = Any[]
    merged = req()

    for (idx, sym) in enumerate(items)
        term = f(sym, idx)
        ae = _to_aexpr(term)
        push!(nodes, ae.node)
        merged = merge_requirements(merged, ae.req)
    end

    return AExpr(sum_expr(nodes), merged)
end

# Defensive swapped-argument-order method for do-block parsing oddities.
@inline sum_over(f::Function, items) = sum_over(items, f)

# -----------------------------------------------------------------------------
# Equation wrapper
# -----------------------------------------------------------------------------

struct Equation
    ae::AExpr
end

@inline expr(eq::Equation) = eq.ae.node
@inline requirements(eq::Equation) = eq.ae.req

end # module Equations
