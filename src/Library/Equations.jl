"""Agate.Library.Equations

Construction-time symbolic building blocks for authoring biogeochemical dynamics
without hand-written `Expr` quoting and without any `@register_*` macros.

This module is intentionally non-circular and does not depend on Models or
Constructor. Objects defined here are used *only at construction time*.

Design notes
------------
- Dynamics authors write equations using *bare parameter identifiers* like
  `maximum_growth_rate[i]` or `detritus_remineralization`.
- Those identifiers are bound (in model modules) to small `ParamVar` placeholders.
- When a `ParamVar` is indexed, it produces an `AExpr` node and records a
  requirement (scalar/vector/matrix) for the parameter key.
- At runtime, parameters are plain numeric scalars and arrays; the symbolic
  placeholders never enter kernels.

Public API:
- `@paramvars a b c` — define `const a = ParamVar{:a}()` placeholders in the
  calling module.
- `declare_parameter_vars!(mod, names; export_vars=true)` — programmatic variant.
- `Σ(items) do sym, idx ... end` — construction-time symbolic sum builder.
- `Equation(::AExpr)` — wrapper returned by dynamics builders.

The old `group_param`/`community_param` distinction is intentionally removed;
missing/`nothing` policy is specified in the parameter registry.
"""

module Equations

using ..ExprUtils: sum_expr
using Logging

export Requirements, req, merge_requirements
export AExpr, Equation, expr, requirements
export ParamVar, @paramvars, declare_parameter_vars!
export Σ

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

"""Define `const` parameter placeholders in the calling module.

Example:

```julia
@paramvars linear_mortality quadratic_mortality

mort = linear_mortality[i] * P
```
"""
macro paramvars(names...)
    defs = Expr[]
    for n in names
        n isa Symbol || error("@paramvars expects bare identifiers")
        push!(defs, :(const $(esc(n)) = ParamVar{$(QuoteNode(n))}()))
    end
    return Expr(:block, defs...)
end

"""Programmatically declare placeholders in `mod` for the given `names`.

If `export_vars=true`, also exports each declared name from `mod`.
"""
function declare_parameter_vars!(mod::Module, names::AbstractVector{Symbol}; export_vars::Bool=true)
    for name in names
        # Skip if already defined to avoid clobbering user bindings.
        if isdefined(mod, name)
            @warn "Parameter var $(name) already defined in $(mod); skipping." maxlog=1
            continue
        end
        Core.eval(mod, :(const $(name) = $(Equations).ParamVar{$(QuoteNode(name))}()))
        if export_vars
            Core.eval(mod, :(export $(name)))
        end
    end
    return nothing
end

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

Base.:+(a::AExpr, b::AExpr) = _binop(:+, a, b)
Base.:-(a::AExpr, b::AExpr) = _binop(:-, a, b)
Base.:*(a::AExpr, b::AExpr) = _binop(:*, a, b)
Base.:/(a::AExpr, b::AExpr) = _binop(:/, a, b)
Base.:^(a::AExpr, b::AExpr) = _binop(:^, a, b)
Base.:-(a::AExpr) = _unop(:-, a)

# -----------------------------------------------------------------------------
# ParamVar arithmetic
#
# ParamVar is used as a symbolic placeholder inside equation DSLs, e.g.
#
#     detritus_remineralization * :D
#     1 - DOM_POM_fractionation
#
# The core expression builder operates on AExpr, but ParamVar is *not* a subtype
# of AExpr. Define a small set of arithmetic overloads so ParamVar can appear on
# either side of standard operators.

const _EquationExprLike = Union{AExpr, ParamVar, Symbol, Number}

Base.:-(a::ParamVar) = _unop(:-, a)

Base.:+(a::ParamVar, b::_EquationExprLike) = _binop(:+, a, b)
Base.:+(a::_EquationExprLike, b::ParamVar) = _binop(:+, a, b)

# Disambiguate ParamVar ⨯ ParamVar (otherwise Julia sees two equally-specific methods).
Base.:+(a::ParamVar, b::ParamVar) = _binop(:+, a, b)

Base.:-(a::ParamVar, b::_EquationExprLike) = _binop(:-, a, b)
Base.:-(a::_EquationExprLike, b::ParamVar) = _binop(:-, a, b)
Base.:-(a::ParamVar, b::ParamVar) = _binop(:-, a, b)

Base.:*(a::ParamVar, b::_EquationExprLike) = _binop(:*, a, b)
Base.:*(a::_EquationExprLike, b::ParamVar) = _binop(:*, a, b)
Base.:*(a::ParamVar, b::ParamVar) = _binop(:*, a, b)

Base.:/(a::ParamVar, b::_EquationExprLike) = _binop(:/, a, b)
Base.:/(a::_EquationExprLike, b::ParamVar) = _binop(:/, a, b)
Base.:/(a::ParamVar, b::ParamVar) = _binop(:/, a, b)

Base.:^(a::ParamVar, b::_EquationExprLike) = _binop(:^, a, b)
Base.:^(a::_EquationExprLike, b::ParamVar) = _binop(:^, a, b)
Base.:^(a::ParamVar, b::ParamVar) = _binop(:^, a, b)



# Disambiguation for ParamVar <-> AExpr binary ops.
# Without these, calls like *(ParamVar, AExpr) are ambiguous between the
# ParamVar overloads above and the generic (a, b::AExpr) fallbacks below.
Base.:+(a::ParamVar, b::AExpr) = _binop(:+, a, b)
Base.:+(a::AExpr, b::ParamVar) = _binop(:+, a, b)

Base.:-(a::ParamVar, b::AExpr) = _binop(:-, a, b)
Base.:-(a::AExpr, b::ParamVar) = _binop(:-, a, b)

Base.:*(a::ParamVar, b::AExpr) = _binop(:*, a, b)
Base.:*(a::AExpr, b::ParamVar) = _binop(:*, a, b)

Base.:/(a::ParamVar, b::AExpr) = _binop(:/, a, b)
Base.:/(a::AExpr, b::ParamVar) = _binop(:/, a, b)

Base.:^(a::ParamVar, b::AExpr) = _binop(:^, a, b)
Base.:^(a::AExpr, b::ParamVar) = _binop(:^, a, b)

# Mixed with literals / symbols / ParamVar
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

# Defensive swapped-argument-order method for do-block parsing oddities.
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
