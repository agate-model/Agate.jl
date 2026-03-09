"""Map tracer, group, and auxiliary names to positional runtime indices.

`TracerIndex` stores only integer offsets and is safe to embed in GPU kernels.
The symbol sets are encoded in type parameters for compile-time specialization.

Type parameters
---------------
- `TR`: tracer symbols in positional order.
- `GS`: plankton group symbols.
- `AF`: auxiliary field symbols.
- `NG`: number of plankton groups.

Fields
------
- `n_tracers`: number of tracer arguments.
- `aux_base`: one-based position of the first auxiliary value in `args`.
- `plankton_base`: one-based position of the first plankton tracer in `args`, or
  `0` when no plankton layout is defined.
- `group_bases`: one-based positions of each group's first class tracer.
- `group_counts`: number of classes in each group.
"""
struct TracerIndex{TR,GS,AF,NG}
    n_tracers::Int
    aux_base::Int
    plankton_base::Int
    group_bases::NTuple{NG,Int}
    group_counts::NTuple{NG,Int}
end

"""Accessor wrapper used inside kernels.

`Tracers` is a lightweight handle that exposes convenient, readable accessors
such as `tracers.P(args, i)` and `tracers.N(args)`.
"""
struct Tracers{TI}
    idx::TI
end

"""Callable accessor returned by `tr.<name>`.

The accessor kind is encoded in the `K` type parameter:
- `:group` for group-level access that takes a class index
- `:scalar` for scalar tracer/aux access with no class index
- `:plankton` for global plankton-class access
"""
struct TracerAccessor{S,TI,K}
    idx::TI
end

@inline (a::TracerAccessor{S,TI,:scalar})(args) where {S,TI} =
    _scalar_value(a.idx, args, Val(S))

@inline (a::TracerAccessor{S,TI,:group})(args, i::Int) where {S,TI} =
    _group_value(a.idx, args, Val(S), i)

@inline (a::TracerAccessor{S,TI,:plankton})(args, i::Int) where {S,TI} =
    _plankton_value(a.idx, args, i)

"""Create a default `TracerIndex` for an arbitrary tracer set.

This constructor does not assume any plankton-group structure.
"""
function build_tracer_index(tracers::Tuple, auxiliary_fields::Tuple)
    TR = tracers
    AF = auxiliary_fields
    return TracerIndex{TR,(),AF,0}(length(tracers), length(tracers) + 1, 0, (), ())
end

"""Create a `TracerIndex` using a parsed community context.

`n_biogeochem_tracers` is the count of non-plankton tracers appearing before
plankton tracers in the positional tracer list.
"""
function build_tracer_index(
    community_context, tracers::Tuple, auxiliary_fields::Tuple; n_biogeochem_tracers::Int
)
    groups_vec = Symbol[]
    if !isempty(community_context.group_symbols)
        last = community_context.group_symbols[1]
        push!(groups_vec, last)
        @inbounds for i in 2:length(community_context.group_symbols)
            g = community_context.group_symbols[i]
            if g !== last
                push!(groups_vec, g)
                last = g
            end
        end
    end

    groups = Tuple(groups_vec)
    NG = length(groups)

    plankton_base = n_biogeochem_tracers + 1
    bases = ntuple(i -> begin
        g = groups[i]
        first_global = first(community_context.group_indices[g])
        plankton_base + first_global - 1
    end, NG)
    counts = ntuple(i -> length(community_context.group_indices[groups[i]]), NG)

    TR = tracers
    AF = auxiliary_fields
    return TracerIndex{TR,groups,AF,NG}(
        length(tracers), length(tracers) + 1, plankton_base, bases, counts
    )
end

@generated function _find_in_tuple(::Val{sym}, ::Val{tup}) where {sym,tup}
    for (i, s) in enumerate(tup)
        if s === sym
            return :($i)
        end
    end
    return :(0)
end

@inline _tracer_position(::TracerIndex{TR}, ::Val{sym}) where {TR,sym} =
    _find_in_tuple(Val(sym), Val(TR))

@inline _aux_slot(::TracerIndex{TR,GS,AF,NG}, ::Val{sym}) where {TR,GS,AF,NG,sym} =
    _find_in_tuple(Val(sym), Val(AF))

@inline _group_slot(::TracerIndex{TR,GS,AF,NG}, ::Val{g}) where {TR,GS,AF,NG,g} =
    _find_in_tuple(Val(g), Val(GS))

@inline function _scalar_value(
    idx::TracerIndex{TR,GS,AF,NG}, args, ::Val{sym}
) where {TR,GS,AF,NG,sym}
    pos = _tracer_position(idx, Val(sym))
    if pos == 0
        slot = _aux_slot(idx, Val(sym))
        slot == 0 && throw(ArgumentError("Unknown tracer/auxiliary name :$sym"))
        return @inbounds args[idx.aux_base + (slot - 1)]
    end
    return @inbounds args[pos]
end

@inline function _group_value(
    idx::TracerIndex{TR,GS,AF,NG}, args, ::Val{g}, i::Int
) where {TR,GS,AF,NG,g}
    slot = _group_slot(idx, Val(g))
    slot == 0 && throw(ArgumentError("Unknown group :$g"))
    base = @inbounds idx.group_bases[slot]
    return @inbounds args[base + (i - 1)]
end

@inline function _plankton_value(idx::TracerIndex, args, i::Int)
    idx.plankton_base == 0 &&
        throw(ArgumentError("No plankton layout is defined for this TracerIndex."))
    return @inbounds args[idx.plankton_base + (i - 1)]
end

@generated function Base.getproperty(tr::Tracers{TI}, name::Symbol) where {TI}
    TR = TI.parameters[1]
    GS = TI.parameters[2]
    AF = TI.parameters[3]

    idx_expr = :(getfield(tr, :idx))
    branches = Tuple{Expr, Expr}[]

    push!(branches, (:(name === :idx), :(return $idx_expr)))
    push!(branches, (:(name === :plankton), :(return TracerAccessor{:plankton, TI, :plankton}($idx_expr))))

    for g in GS
        g_q = QuoteNode(g)
        push!(branches, (:(name === $g_q), :(return TracerAccessor{$g_q, TI, :group}($idx_expr))))
    end

    for s in TR
        s_q = QuoteNode(s)
        push!(branches, (:(name === $s_q), :(return TracerAccessor{$s_q, TI, :scalar}($idx_expr))))
    end

    for a in AF
        a_q = QuoteNode(a)
        push!(branches, (:(name === $a_q), :(return TracerAccessor{$a_q, TI, :scalar}($idx_expr))))
    end
    ex = :(throw(ArgumentError("Unknown tracer/group/auxiliary name")))
    for i in reverse(eachindex(branches))
        cond, body = branches[i]
        ex = Expr(:if, cond, body, ex)
    end

    return ex
end
