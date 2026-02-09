import Adapt

"""Rectangular consumer-by-prey interaction matrices and axis mappings.

`InteractionMatrices` is the canonical runtime representation for role-aware
interaction matrices.

Interaction data is stored in rectangular matrices sized `(n_consumer, n_prey)`
where:

- `n_consumer = length(community_context.consumer_indices)`
- `n_prey     = length(community_context.prey_indices)`

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
    community_context::CommunityContext,
    value,
    row_indices::Vector{Int},
    col_indices::Vector{Int},
    key::Symbol,
)
    n_total = community_context.n_total
    nr = length(row_indices)
    nc = length(col_indices)
    if value isa AbstractMatrix
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
                "interaction matrix '$key' must be a matrix; got $(typeof(value))"
            ),
        )
    end
end

function finalize_interaction_parameters(
    factory::AbstractBGCFactory, community_context::CommunityContext, params::NamedTuple
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

    consumer_indices = community_context.consumer_indices
    prey_indices = community_context.prey_indices

    pal_rect = _rect_value_for_axes(
        community_context,
        params.palatability_matrix,
        consumer_indices,
        prey_indices,
        :palatability_matrix,
    )
    assim_rect = _rect_value_for_axes(
        community_context,
        params.assimilation_matrix,
        consumer_indices,
        prey_indices,
        :assimilation_matrix,
    )

    global_to_consumer = _inverse_axis_map(consumer_indices, community_context.n_total)
    global_to_prey = _inverse_axis_map(prey_indices, community_context.n_total)

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

"""Normalize `interaction_overrides` into a `NamedTuple` of parameter overrides.

`interaction_overrides` may be:

- `nothing` (no overrides)
- a `NamedTuple` of updates

Interaction overrides are **data-only**. Values must be explicit, canonical
axis-sized rectangular matrices. Provider functions / callables are **not**
supported.

For a matrix parameter with declared `axes` (for example `(:consumer, :prey)`),
users must pass a rectangular matrix sized to the declared axes (for example
`(n_consumer, n_prey)`).

If you need to derive matrices from traits or other parameters, define a
`Variant` / `Factory` default that produces concrete rectangular matrices during
construction.
"""
function normalize_interaction_overrides(
    factory::AbstractBGCFactory,
    community_context::CommunityContext{FT},
    interaction_overrides::Union{Nothing,NamedTuple},
) where {FT}
    interaction_overrides === nothing && return (;)

    resolved = ()
    for (key, value) in pairs(interaction_overrides)
        spec = parameter_spec(factory, key)
        spec !== nothing || throw(
            ArgumentError(
                "interaction override '$key' is missing a ParameterSpec in parameter_directory(::$(typeof(factory))).",
            ),
        )

        (spec.shape === :matrix) || throw(
            ArgumentError(
                "interaction override '$key' must target a matrix ParameterSpec; got shape $(spec.shape). Use `parameters=(; ...)` for non-matrix parameters.",
            ),
        )
        (spec.axes !== nothing) || throw(
            ArgumentError(
                "interaction override '$key' must target an axes-declared matrix ParameterSpec (canonical rectangular matrices only). Providers are not supported; use a Variant/Factory if you need derived matrices.",
            ),
        )

        row_axis, col_axis = spec.axes
        row_indices = axis_indices(community_context, row_axis)
        col_indices = axis_indices(community_context, col_axis)
        nr = length(row_indices)
        nc = length(col_indices)

        if value isa Function || applicable(value, community_context)
            throw(
                ArgumentError(
                    "interaction override '$key' expected rectangular matrix of size ($(nr), $(nc)) for axes $(spec.axes); providers are not supported; use a Variant/Factory if you need derived matrices.",
                ),
            )
        end

        value isa AbstractMatrix || throw(
            ArgumentError(
                "interaction override '$key' expected rectangular matrix of size ($(nr), $(nc)) for axes $(spec.axes); got $(typeof(value)). Providers are not supported; use a Variant/Factory if you need derived matrices.",
            ),
        )

        resolved_value = _normalize_interaction_value(community_context, spec, key, value)
        _validate_interaction_override(community_context, spec, key, resolved_value)
        resolved = (resolved..., key => resolved_value)
    end

    return (; resolved...)
end


"""Return the global plankton indices for an interaction axis.

Axes may be:

- `:consumer` (role-defined consumer axis)
- `:prey` (role-defined prey axis)
- any existing group `Symbol` present in `community_context.group_indices`
"""
@inline axis_indices(community_context::CommunityContext, axis::Symbol) =
    if axis === :consumer
        community_context.consumer_indices
    elseif axis === :prey
        community_context.prey_indices
    elseif haskey(community_context.group_indices, axis)
        community_context.group_indices[axis]
    else
        throw(
            ArgumentError(
                "Unknown interaction axis '$axis'. Valid axes are :consumer, :prey, or an existing group symbol.",
            ),
        )
    end

@inline function _normalize_interaction_value(
    community_context::CommunityContext, spec::ParameterSpec, key::Symbol, value
)
    spec.shape === :matrix || return value
    n_total = community_context.n_total

    groups = unique(community_context.group_symbols)
    ng = length(groups)

    # Axes-aware matrices are normalized into their axis-local rectangular form.
    if spec.axes !== nothing
        row_axis, col_axis = spec.axes
        row_indices = axis_indices(community_context, row_axis)
        col_indices = axis_indices(community_context, col_axis)

        nr = length(row_indices)
        nc = length(col_indices)

        value isa AbstractMatrix || return value

        if size(value) == (nr, nc)
            return _convert_matrix_eltype(community_context, value)
        end

        if size(value) == (ng, ng)
            throw(
                ArgumentError(
                    "interaction override '$key' looks like a group-by-group matrix (size $(ng)x$(ng)). Group-block overrides are not supported; pass a $(nr)x$(nc) axis matrix for axes $(spec.axes). Providers are not supported; use a Variant/Factory if you need derived matrices.",
                ),
            )
        end

        throw(
            ArgumentError(
                "interaction override '$key' must be a $(nr)x$(nc) axis matrix for axes $(spec.axes); got size $(size(value))",
            ),
        )
    end
    # Non-axes matrices normalize to full square matrices.

    value isa AbstractMatrix || return value

    size(value) == (n_total, n_total) && return _convert_matrix_eltype(community_context, value)

    if size(value) == (ng, ng)
        throw(
            ArgumentError(
                "interaction override '$key' looks like a group-by-group matrix (size $(ng)x$(ng)). Group-block overrides are not supported; pass a $(n_total)x$(n_total) matrix instead.",
            ),
        )
    end

    throw(
        ArgumentError(
            "interaction override '$key' must be a $(n_total)x$(n_total) matrix; got size $(size(value))",
        ),
    )
end


@inline function _convert_matrix_eltype(community_context::CommunityContext, A::AbstractMatrix)
    eltype(A) === community_context.FT && return A
    R = similar(A, community_context.FT, size(A)...)
    copyto!(R, A)
    return R
end

@inline function _validate_interaction_override(
    community_context::CommunityContext, spec::ParameterSpec, key::Symbol, value
)
    if spec.shape === :matrix
        value isa AbstractMatrix || throw(
            ArgumentError(
                "interaction override '$key' must be a matrix; got $(typeof(value))"
            ),
        )

        if spec.axes === nothing
            n_total = community_context.n_total
            size(value) == (n_total, n_total) || throw(
                ArgumentError(
                    "interaction override '$key' must be a $(n_total)x$(n_total) matrix after normalization; got size $(size(value))",
                ),
            )
        else
            row_axis, col_axis = spec.axes
            row_indices = axis_indices(community_context, row_axis)
            col_indices = axis_indices(community_context, col_axis)
            nr = length(row_indices)
            nc = length(col_indices)
            size(value) == (nr, nc) || throw(
                ArgumentError(
                    "interaction override '$key' must be a $(nr)x$(nc) matrix for axes $(spec.axes) after normalization; got size $(size(value))",
                ),
            )
        end

    elseif spec.shape === :vector
        value isa AbstractVector || throw(
            ArgumentError(
                "interaction override '$key' must be a vector; got $(typeof(value))"
            ),
        )

        n_total = community_context.n_total
        length(value) == n_total || throw(
            ArgumentError(
                "interaction override '$key' must have length $n_total (got $(length(value)))",
            ),
        )
    end

    return nothing
end
