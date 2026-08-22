"""Return parameter specs stored on the canonical consumer-by-prey axes."""
function interaction_parameter_specs(source)
    return Tuple(
        spec for spec in parameter_directory(source) if
        spec.shape === :matrix && spec.axes == (:consumer, :prey)
    )
end

"""Return runtime metadata for consumer-by-prey parameter matrices.

Interaction matrices themselves remain ordinary top-level parameters. This metadata
stores only the canonical parameter names and the ecological class labels used by the
consumer and prey axes, so introspection and active-parameter selection do not require
a second parameter representation.
"""
function interaction_axis_metadata(source, community_context::CommunityContext)
    specs = interaction_parameter_specs(source)
    isempty(specs) && return nothing

    consumer_classes = Tuple(
        community_context.class_symbols[i] for i in community_context.consumer_indices
    )
    prey_classes = Tuple(
        community_context.class_symbols[i] for i in community_context.prey_indices
    )
    parameter_names = Tuple(spec.name for spec in specs)

    return (;
        parameters=parameter_names,
        consumers=consumer_classes,
        prey=prey_classes,
    )
end

"""Return the global plankton indices for a semantic storage axis.

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
