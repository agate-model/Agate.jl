module Functors

export Requirements, req, merge_requirements, CompiledEquation, requirements

"""A small immutable set of parameter requirements.

The requirement lists are tuples of `Symbol`s.
"""
struct Requirements{V,M,S}
    vectors::V
    matrices::M
    scalars::S
end

_symbol_tuple(::Nothing) = ()
_symbol_tuple(x::Symbol) = (x,)
_symbol_tuple(x::Tuple) = Symbol.(x)
_symbol_tuple(x::AbstractVector) = Tuple(Symbol.(x))
function _symbol_tuple(x)
    throw(
        ArgumentError(
            "Requirements must be Symbol, Tuple, Vector, or nothing; got $(typeof(x))."
        ),
    )
end

"""    req(; vectors=(), matrices=(), scalars=())

Construct a `Requirements` object.

Each keyword accepts a `Symbol`, a tuple of symbols, a vector of symbols, or `nothing`.
"""
function req(; vectors=(), matrices=(), scalars=())
    return Requirements(
        _symbol_tuple(vectors), _symbol_tuple(matrices), _symbol_tuple(scalars)
    )
end

"""A wrapper that associates a callable with its parameter requirements."""
struct CompiledEquation{F,R}
    f::F
    r::R
end

"""Return the parameter requirements for a compiled equation."""
@inline requirements(eq::CompiledEquation) = eq.r

function _merge_syms(a::Tuple, b::Tuple)
    isempty(b) && return a
    isempty(a) && return b

    out = Symbol[]
    for s in a
        push!(out, s)
    end
    for s in b
        (s in out) || push!(out, s)
    end

    return Tuple(out)
end

"""    merge_requirements(r1, r2)

Combine two `Requirements` objects.
"""
function merge_requirements(r1::Requirements, r2::Requirements)
    return req(;
        vectors=_merge_syms(r1.vectors, r2.vectors),
        matrices=_merge_syms(r1.matrices, r2.matrices),
        scalars=_merge_syms(r1.scalars, r2.scalars),
    )
end

end # module
