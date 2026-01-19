import Adapt

"""Rectangular consumer-by-prey interaction matrices and axis mappings.

`InteractionMatrices` is the canonical runtime representation for role-aware
interaction matrices. For compatibility with the existing tracer code, the
constructor pipeline also builds `InteractionMatrixView` objects that present
the rectangular matrices as zero-padded square `(n_total, n_total)` matrices.
"""

struct InteractionMatrices{PM, AM, VI, MI}
    palatability::PM
    assimilation::AM
    consumer_global::VI
    prey_global::VI
    global_to_consumer::MI
    global_to_prey::MI
end
Adapt.@adapt_structure InteractionMatrices

"""A zero-padded square view of a rectangular interaction matrix.

This wrapper supports the legacy access pattern `M[predator_idx, prey_idx]` in
tracer kernels while storing the interaction data in a rectangular
consumer-by-prey matrix.
"""
struct InteractionMatrixView{T, M, MI} <: AbstractMatrix{T}
    rect::M
    global_to_row::MI
    global_to_col::MI
end

InteractionMatrixView(rect::M, global_to_row::MI, global_to_col::MI) where {M, MI} =
    InteractionMatrixView{eltype(rect), M, MI}(rect, global_to_row, global_to_col)

Adapt.@adapt_structure InteractionMatrixView

Base.eltype(V::InteractionMatrixView) = eltype(V.rect)
Base.size(V::InteractionMatrixView) = (length(V.global_to_row), length(V.global_to_col))
Base.IndexStyle(::Type{<:InteractionMatrixView}) = IndexCartesian()

@inline function Base.getindex(V::InteractionMatrixView, i::Int, j::Int)
    ri = @inbounds V.global_to_row[i]
    cj = @inbounds V.global_to_col[j]
    if ri == 0 || cj == 0
        return zero(eltype(V))
    end
    return @inbounds V.rect[ri, cj]
end

@inline function _inverse_axis_map(axis_indices::AbstractVector{Int}, n_total::Int)
    m = zeros(Int, n_total)
    # NOTE: avoid reserved keywords (`global`, `local`) as variable names
    for (lidx, gidx) in pairs(axis_indices)
        @inbounds m[gidx] = lidx
    end
    return m
end

@inline function _rect_value_for_axes(
    ctx::InteractionContext,
    value,
    row_indices::Vector{Int},
    col_indices::Vector{Int},
    key::Symbol,
)
    n_total = ctx.n_total
    nr = length(row_indices)
    nc = length(col_indices)

    if value isa InteractionMatrixView
        return value.rect
    elseif value isa GroupBlockMatrix
        full = expand_group_block_matrix(ctx, value.B)
        return full[row_indices, col_indices]
    elseif value isa AbstractMatrix
        if size(value) == (nr, nc)
            return value
        elseif size(value) == (n_total, n_total)
            return value[row_indices, col_indices]
        else
            throw(ArgumentError(
                "interaction matrix '$key' must be $(nr)x$(nc) (axes) or $(n_total)x$(n_total) (full); got size $(size(value))",
            ))
        end
    else
        throw(ArgumentError("interaction matrix '$key' must be a matrix; got $(typeof(value))"))
    end
end

function finalize_interaction_parameters(factory::AbstractBGCFactory, ctx::InteractionContext, params::NamedTuple)
    spec_pal = parameter_spec(factory, :palatability_matrix)
    spec_assim = parameter_spec(factory, :assimilation_matrix)

    if spec_pal === nothing || spec_assim === nothing
        return params
    end
    if spec_pal.axes != (:consumer, :prey) || spec_assim.axes != (:consumer, :prey)
        return params
    end
    if !haskey(params, :palatability_matrix) || !haskey(params, :assimilation_matrix)
        return params
    end

    consumer_indices = ctx.consumer_indices
    prey_indices = ctx.prey_indices

    pal_rect = _rect_value_for_axes(ctx, params.palatability_matrix, consumer_indices, prey_indices, :palatability_matrix)
    assim_rect = _rect_value_for_axes(ctx, params.assimilation_matrix, consumer_indices, prey_indices, :assimilation_matrix)

    global_to_consumer = _inverse_axis_map(consumer_indices, ctx.n_total)
    global_to_prey = _inverse_axis_map(prey_indices, ctx.n_total)

    pal_view = InteractionMatrixView(pal_rect, global_to_consumer, global_to_prey)
    assim_view = InteractionMatrixView(assim_rect, global_to_consumer, global_to_prey)

    interactions = InteractionMatrices(
        pal_rect,
        assim_rect,
        consumer_indices,
        prey_indices,
        global_to_consumer,
        global_to_prey,
    )

    return merge(
        params,
        (palatability_matrix=pal_view, assimilation_matrix=assim_view, interactions=interactions),
    )
end
