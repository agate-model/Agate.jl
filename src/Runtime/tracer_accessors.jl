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

@inline (a::TracerAccessor{S,TI,:scalar})(args) where {S,TI} = scalar_value(
    a.idx, args, Val(S)
)

@inline (a::TracerAccessor{S,TI,:group})(args, i::Int) where {S,TI} = group_value(
    a.idx, args, Val(S), i
)

@inline (a::TracerAccessor{S,TI,:plankton})(args, i::Int) where {S,TI} = plankton_value(
    a.idx, args, i
)

"""Create a default `TracerIndex` for an arbitrary tracer set.

This constructor does not assume any plankton-group structure.
"""
function build_tracer_index(tracers::Tuple, auxiliary_fields::Tuple)
    TR = tracers
    AF = auxiliary_fields
    return TracerIndex{TR,(),AF,0}(length(tracers), length(tracers) + 1, 0, (), ())
end

"""Create a `TracerIndex` directly from one realized `ModelLayout`.

Single-state population classes expose the familiar group and global-plankton accessors.
For multi-state populations ecological classes do not map one-to-one onto physical tracer
arguments, so state-aware compiled processes use their pre-resolved tracer operands and the
runtime accessor remains scalar-only.
"""
function build_tracer_index(layout::ModelLayout)
    tracers = layout.tracer_order
    auxiliary_fields = layout.auxiliary_fields
    class_positions = Tuple(
        getproperty(layout.tracer_indices, class)
        for class in layout.class_symbols if hasproperty(layout.tracer_indices, class)
    )
    if length(class_positions) != length(layout.class_symbols)
        return build_tracer_index(tracers, auxiliary_fields)
    end

    groups = keys(layout.group_indices)
    NG = length(groups)
    bases = ntuple(NG) do i
        first_class = layout.class_symbols[first(getproperty(layout.group_indices, groups[i]))]
        getproperty(layout.tracer_indices, first_class)
    end
    counts = ntuple(i -> length(getproperty(layout.group_indices, groups[i])), NG)
    plankton_base = isempty(class_positions) ? 0 : first(class_positions)

    return TracerIndex{tracers,groups,auxiliary_fields,NG}(
        length(tracers), length(tracers) + 1, plankton_base, bases, counts
    )
end

@generated function find_in_tuple(::Val{sym}, ::Val{tup}) where {sym,tup}
    for (i, s) in enumerate(tup)
        if s === sym
            return :($i)
        end
    end
    return :(0)
end

@inline tracer_position(::TracerIndex{TR}, ::Val{sym}) where {TR,sym} = find_in_tuple(
    Val(sym), Val(TR)
)

@inline aux_slot(::TracerIndex{TR,GS,AF,NG}, ::Val{sym}) where {TR,GS,AF,NG,sym} = find_in_tuple(
    Val(sym), Val(AF)
)

@inline group_slot(::TracerIndex{TR,GS,AF,NG}, ::Val{g}) where {TR,GS,AF,NG,g} = find_in_tuple(
    Val(g), Val(GS)
)

@inline function scalar_value(
    idx::TracerIndex{TR,GS,AF,NG}, args, ::Val{sym}
) where {TR,GS,AF,NG,sym}
    pos = tracer_position(idx, Val(sym))
    if pos == 0
        slot = aux_slot(idx, Val(sym))
        slot == 0 && throw(ArgumentError("Unknown tracer/auxiliary name :$sym"))
        return @inbounds args[idx.aux_base + (slot - 1)]
    end
    return @inbounds args[pos]
end

@inline function group_value(
    idx::TracerIndex{TR,GS,AF,NG}, args, ::Val{g}, i::Int
) where {TR,GS,AF,NG,g}
    slot = group_slot(idx, Val(g))
    slot == 0 && throw(ArgumentError("Unknown group :$g"))
    base = @inbounds idx.group_bases[slot]
    return @inbounds args[base + (i - 1)]
end

@inline function plankton_value(idx::TracerIndex, args, i::Int)
    idx.plankton_base == 0 &&
        throw(ArgumentError("No plankton layout is defined for this TracerIndex."))
    return @inbounds args[idx.plankton_base + (i - 1)]
end

@generated function Base.getproperty(tr::Tracers{TI}, name::Symbol) where {TI}
    TR = TI.parameters[1]
    GS = TI.parameters[2]
    AF = TI.parameters[3]

    idx_expr = :(getfield(tr, :idx))
    branches = Tuple{Expr,Expr}[]

    push!(branches, (:(name === :idx), :(return $idx_expr)))
    push!(
        branches,
        (
            :(name === :plankton),
            :(return TracerAccessor{:plankton,TI,:plankton}($idx_expr)),
        ),
    )

    for g in GS
        g_q = QuoteNode(g)
        push!(
            branches,
            (:(name === $g_q), :(return TracerAccessor{$g_q,TI,:group}($idx_expr))),
        )
    end

    for s in TR
        s_q = QuoteNode(s)
        push!(
            branches,
            (:(name === $s_q), :(return TracerAccessor{$s_q,TI,:scalar}($idx_expr))),
        )
    end

    for a in AF
        a_q = QuoteNode(a)
        push!(
            branches,
            (:(name === $a_q), :(return TracerAccessor{$a_q,TI,:scalar}($idx_expr))),
        )
    end
    ex = :(throw(ArgumentError("Unknown tracer/group/auxiliary name")))
    for i in reverse(eachindex(branches))
        cond, body = branches[i]
        ex = Expr(:if, cond, body, ex)
    end

    return ex
end
