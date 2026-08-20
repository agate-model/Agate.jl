import Adapt

"""Rectangular consumer-by-prey interaction matrices and axis mappings.

`InteractionMatrices` is the canonical runtime representation for role-aware
interaction matrices. Arbitrary consumer-by-prey matrices are stored in the
concrete `NamedTuple` field `matrices`; the remaining fields describe the
shared consumer and prey axes.

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

Matrix keys use each parameter spec's explicit runtime path. A consumer-by-prey
parameter stored at `(:interactions, :encounter)` is available as
`interactions.encounter`.
"""
struct InteractionMatrices{M,VI1,VI2,MI1,MI2}
    matrices::M
    consumer_global::VI1
    prey_global::VI2
    global_to_consumer::MI1
    global_to_prey::MI2
end
Adapt.@adapt_structure InteractionMatrices

@inline function Base.getproperty(interactions::InteractionMatrices, name::Symbol)
    if name === :matrices ||
       name === :consumer_global ||
       name === :prey_global ||
       name === :global_to_consumer ||
       name === :global_to_prey
        return getfield(interactions, name)
    end

    matrices = getfield(interactions, :matrices)
    hasproperty(matrices, name) && return getproperty(matrices, name)
    return getfield(interactions, name)
end

@inline function Base.propertynames(interactions::InteractionMatrices, private::Bool=false)
    structural = fieldnames(typeof(interactions))
    matrix_names = propertynames(getfield(interactions, :matrices), private)
    return (structural..., matrix_names...)
end

"""Compare resolved interaction state structurally for deterministic replay.

Every stored field participates in equality. `isequal` preserves exact
comparison for matching non-finite values and for arbitrary named interaction
matrices.
"""
function Base.:(==)(a::InteractionMatrices, b::InteractionMatrices)
    fieldcount(typeof(a)) == fieldcount(typeof(b)) || return false
    return all(i -> isequal(getfield(a, i), getfield(b, i)), 1:fieldcount(typeof(a)))
end

@inline function inverse_axis_map(axis_indices::AbstractVector{Int}, n_total::Int)
    m = zeros(Int, n_total)
    for (lidx, gidx) in pairs(axis_indices)
        @inbounds m[gidx] = lidx
    end
    return m
end

@inline function rect_value_for_axes(
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

function interaction_parameter_specs(source)
    return Tuple(
        spec for spec in parameter_directory(source) if
        spec.shape === :matrix && spec.axes == (:consumer, :prey)
    )
end

interaction_parameter_names(source) =
    Tuple(spec.name for spec in interaction_parameter_specs(source))

function interaction_runtime_name(spec::ParameterSpec)
    path = spec.runtime_path
    length(path) == 2 && first(path) === :interactions || throw(
        ArgumentError(
            "consumer-by-prey parameter :$(spec.name) must declare runtime_path=(:interactions, :name)",
        ),
    )
    return last(path)
end

function finalize_interaction_parameters(
    source, community_context::CommunityContext, params::NamedTuple
)
    specs = interaction_parameter_specs(source)
    isempty(specs) && return params
    parameter_names = Tuple(spec.name for spec in specs)
    all(name -> haskey(params, name), parameter_names) || return params

    runtime_names = Tuple(interaction_runtime_name(spec) for spec in specs)
    length(unique(runtime_names)) == length(runtime_names) || throw(
        ArgumentError(
            "consumer-by-prey interaction parameter names map to duplicate runtime names: $(runtime_names)",
        ),
    )

    consumer_indices = community_context.consumer_indices
    prey_indices = community_context.prey_indices

    rectangular_values = Tuple(
        rect_value_for_axes(
            community_context,
            getproperty(params, name),
            consumer_indices,
            prey_indices,
            name,
        ) for name in parameter_names
    )
    rectangular_parameters = NamedTuple{parameter_names}(rectangular_values)
    matrices = NamedTuple{runtime_names}(rectangular_values)

    global_to_consumer = inverse_axis_map(consumer_indices, community_context.n_total)
    global_to_prey = inverse_axis_map(prey_indices, community_context.n_total)

    interactions = InteractionMatrices(
        matrices,
        consumer_indices,
        prey_indices,
        global_to_consumer,
        global_to_prey,
    )

    return merge(params, rectangular_parameters, (; interactions))
end

"""Normalize `interaction_overrides` into a `NamedTuple` of parameter overrides.

`interaction_overrides` is a `NamedTuple` of updates. Use an empty named tuple
when there are no explicit interaction overrides.

Interaction overrides are **data-only**. Values must be explicit, canonical
axis-sized rectangular matrices.

For a matrix parameter with declared `axes` (for example `(:consumer, :prey)`),
users must pass a rectangular matrix sized to the declared axes (for example
`(n_consumer, n_prey)`).

If you need to derive matrices from traits or other parameters, declare a
`DerivedDefault` in the parameter definition so construction materializes the
concrete rectangular matrix.
"""
function normalize_interaction_overrides(
    source,
    community_context::CommunityContext{T},
    interaction_overrides::NamedTuple,
) where {T}
    resolved = ()
    for (key, value) in pairs(interaction_overrides)
        spec = parameter_spec(source, key)
        spec !== nothing || throw(
            ArgumentError(
                "interaction override '$key' is missing a ParameterSpec in parameter_directory(::$(typeof(source))).",
            ),
        )

        (spec.shape === :matrix) || throw(
            ArgumentError(
                "interaction override '$key' must target a matrix ParameterSpec; got shape $(spec.shape). Use `parameters=(; ...)` for non-matrix parameters.",
            ),
        )
        (spec.axes !== nothing) || throw(
            ArgumentError(
                "interaction override '$key' must target an axes-declared matrix ParameterSpec (canonical rectangular matrices only). Providers are not supported here; declare a DerivedDefault in parameter_definitions for derived values.",
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
                    "interaction override '$key' expected rectangular matrix of size ($(nr), $(nc)) for axes $(spec.axes); providers are not supported here; declare a DerivedDefault in parameter_definitions for derived values.",
                ),
            )
        end

        value isa AbstractMatrix || throw(
            ArgumentError(
                "interaction override '$key' expected rectangular matrix of size ($(nr), $(nc)) for axes $(spec.axes); got $(typeof(value)). Providers are not supported here; declare a DerivedDefault in parameter_definitions for derived values.",
            ),
        )

        resolved_value = normalize_interaction_value(community_context, spec, key, value)
        validate_interaction_override(community_context, spec, key, resolved_value)
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

@inline function normalize_interaction_value(
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
            return convert_matrix_eltype(community_context, value)
        end

        if size(value) == (ng, ng)
            throw(
                ArgumentError(
                    "interaction override '$key' looks like a group-by-group matrix (size $(ng)x$(ng)). Group-block overrides are not supported; pass a $(nr)x$(nc) axis matrix for axes $(spec.axes). Providers are not supported here; declare a DerivedDefault in parameter_definitions for derived values.",
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

    size(value) == (n_total, n_total) &&
        return convert_matrix_eltype(community_context, value)

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

@inline function convert_matrix_eltype(
    community_context::CommunityContext, A::AbstractMatrix
)
    eltype(A) === community_context.scalar_type && return A
    R = similar(A, community_context.scalar_type, size(A)...)
    copyto!(R, A)
    return R
end

@inline function validate_interaction_override(
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
