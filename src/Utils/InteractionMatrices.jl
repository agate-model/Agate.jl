import Adapt

"""Rectangular consumer-by-prey interaction matrices and axis mappings.

`InteractionMatrices` is the canonical runtime representation for role-aware
interaction matrices.

Interaction data is stored in rectangular matrices sized `(n_consumer, n_prey)`
where:

- `n_consumer = length(ctx.consumer_indices)`
- `n_prey     = length(ctx.prey_indices)`

The axis vectors map axis-local indices to global plankton class indices:

- `consumer_global[ic]` gives the global index for consumer axis position `ic`
- `prey_global[ip]` gives the global index for prey axis position `ip`

The inverse maps support fast lookup of axis-local indices from global indices:

- `global_to_consumer[g]` returns `ic` or `0` if `g` is not a consumer
- `global_to_prey[g]` returns `ip` or `0` if `g` is not a prey

No square matrices (or square views) are created or stored.
"""

struct InteractionMatrices{PM,AM,VI,MI}
    palatability::PM
    assimilation::AM
    consumer_global::VI
    prey_global::VI
    global_to_consumer::MI
    global_to_prey::MI
end
Adapt.@adapt_structure InteractionMatrices

@inline function _inverse_axis_map(axis_indices::AbstractVector{Int}, n_total::Int)
    m = zeros(Int, n_total)
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

    if value isa GroupBlockMatrix
        full = expand_group_block_matrix(ctx, value.B)
        return full[row_indices, col_indices]
    elseif value isa AbstractMatrix
        if size(value) == (nr, nc)
            return value
        elseif size(value) == (n_total, n_total)
            return value[row_indices, col_indices]
        else
            throw(
                ArgumentError(
                    "interaction matrix '$key' must be $(nr)x$(nc) (axes) or $(n_total)x$(n_total) (full); got size $(size(value))",
                ),
            )
        end
    else
        throw(
            ArgumentError(
                "interaction matrix '$key' must be a matrix; got $(typeof(value))",
            ),
        )
    end
end

function finalize_interaction_parameters(
    factory::AbstractBGCFactory, ctx::InteractionContext, params::NamedTuple
)
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

    pal_rect = _rect_value_for_axes(
        ctx,
        params.palatability_matrix,
        consumer_indices,
        prey_indices,
        :palatability_matrix,
    )
    assim_rect = _rect_value_for_axes(
        ctx,
        params.assimilation_matrix,
        consumer_indices,
        prey_indices,
        :assimilation_matrix,
    )

    global_to_consumer = _inverse_axis_map(consumer_indices, ctx.n_total)
    global_to_prey = _inverse_axis_map(prey_indices, ctx.n_total)

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
        (
            palatability_matrix=pal_rect,
            assimilation_matrix=assim_rect,
            interactions=interactions,
        ),
    )
end
