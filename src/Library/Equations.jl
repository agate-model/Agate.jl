"""Agate.Library.Equations

Construction-time symbolic building blocks for authoring biogeochemical dynamics
without hand-written `Expr` quoting and without any `@register_*` macros.

This module is intentionally non-circular and does not depend on Models or
Constructor. Objects defined here are used *only at construction time*.

Set C public API (explicit; no aliases):
- group_param(:key)
- community_param(:key)
- interaction_matrix(:key)
- bgc_param(:key)
- Σ(items) do sym, idx ... end
- Equation(::AExpr)
"""

module Equations

using ..ExprUtils: sum_expr

export Requirements, req, merge_requirements
export AExpr, Equation, expr, requirements
export group_param, community_param, interaction_matrix, bgc_param
export Σ

# -----------------------------------------------------------------------------
# Requirements
# -----------------------------------------------------------------------------

"""Construction-time dependency set for a symbolic expression."""
struct Requirements
    group_params::Vector{Symbol}
    community_params::Vector{Symbol}
    matrices::Vector{Symbol}
    scalars::Vector{Symbol}
end

"""Create a `Requirements` object."""
function req(; group=(), community=(), matrices=(), scalars=())
    return Requirements(
        Symbol[group...],
        Symbol[community...],
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
        _append_unique!(copy(r1.group_params), r2.group_params),
        _append_unique!(copy(r1.community_params), r2.community_params),
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

@inline function _to_aexpr(x)
    return x isa AExpr ? x : AExpr(x, req())
end

# -----------------------------------------------------------------------------
# References: parameters, matrices, scalars
# -----------------------------------------------------------------------------

struct ParamRef
    kind::Symbol   # :group or :community
    key::Symbol
end

struct MatrixRef
    key::Symbol
end

"""Group-owned parameter reference (missing in that group => error at construct time)."""
@inline group_param(key::Symbol) = ParamRef(:group, key)

"""Community-optional parameter reference (missing/nothing => 0 at construct time)."""
@inline community_param(key::Symbol) = ParamRef(:community, key)

"""Interaction matrix reference."""
@inline interaction_matrix(key::Symbol) = MatrixRef(key)

"""BGC scalar parameter reference (from `BiogeochemistrySpecification`)."""
@inline function bgc_param(key::Symbol)
    return AExpr(key, req(scalars=(key,)))
end

function Base.getindex(r::ParamRef, i)
    if r.kind === :group
        return AExpr(Expr(:ref, r.key, i), req(group=(r.key,)))
    else
        return AExpr(Expr(:ref, r.key, i), req(community=(r.key,)))
    end
end

function Base.getindex(M::MatrixRef, j, i)
    return AExpr(Expr(:ref, M.key, j, i), req(matrices=(M.key,)))
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

Base.:+(a::AExpr, b::AExpr) = _binop(:+, a, b)
Base.:-(a::AExpr, b::AExpr) = _binop(:-, a, b)
Base.:*(a::AExpr, b::AExpr) = _binop(:*, a, b)
Base.:/(a::AExpr, b::AExpr) = _binop(:/, a, b)
Base.:^(a::AExpr, b::AExpr) = _binop(:^, a, b)
Base.:-(a::AExpr) = _unop(:-, a)

# Mixed with literals / symbols
Base.:+(a::AExpr, b) = _binop(:+, a, b)
Base.:+(a, b::AExpr) = _binop(:+, a, b)
Base.:-(a::AExpr, b) = _binop(:-, a, b)
Base.:-(a, b::AExpr) = _binop(:-, a, b)
Base.:*(a::AExpr, b) = _binop(:*, a, b)
Base.:*(a, b::AExpr) = _binop(:*, a, b)
Base.:/(a::AExpr, b) = _binop(:/, a, b)
Base.:/(a, b::AExpr) = _binop(:/, a, b)
Base.:^(a::AExpr, b) = _binop(:^, a, b)
Base.:^(a, b::AExpr) = _binop(:^, a, b)

# -----------------------------------------------------------------------------
# Σ sum builder (construction-time only)
# -----------------------------------------------------------------------------

"""\
    Σ(items) do sym, idx
        ... -> AExpr
    end

Build a symbolic unrolled sum over `items` (construction-time only).
"""
function Σ(items, f::Function)
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

# Some Julia parsing patterns can end up calling Σ(f, items) when a `do` block is
# attached in unexpected ways. Accept the swapped argument order defensively.
@inline Σ(f::Function, items) = Σ(items, f)

# -----------------------------------------------------------------------------
# Equation wrapper
# -----------------------------------------------------------------------------

struct Equation
    ae::AExpr
end

@inline expr(eq::Equation) = eq.ae.node
@inline requirements(eq::Equation) = eq.ae.req

end # module Equations
