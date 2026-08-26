"""Return runtime metadata for parameters stored on canonical consumer-by-prey axes."""
function interaction_axis_metadata(plan)
    names = Tuple(
        name for name in keys(plan.parameters)
        if begin
            parameter = getproperty(plan.parameters, name)
            parameter.runtime_bound && parameter.storage_axes == (:consumer, :prey)
        end
    )
    isempty(names) && return nothing

    first_parameter = getproperty(plan.parameters, first(names))
    consumers, prey = first_parameter.storage_labels
    for name in Base.tail(names)
        parameter = getproperty(plan.parameters, name)
        parameter.storage_labels == (consumers, prey) || throw(ArgumentError(
            "consumer-by-prey parameters do not share one realized interaction axis",
        ))
    end
    return (; parameters=names, consumers, prey)
end

"""Return global ecological-class indices for one semantic storage axis."""
@inline axis_indices(layout::ModelLayout, axis::Symbol) =
    if axis === :consumer
        layout.consumer_indices
    elseif axis === :prey
        layout.prey_indices
    elseif hasproperty(layout.group_indices, axis)
        getproperty(layout.group_indices, axis)
    else
        throw(
            ArgumentError(
                "Unknown interaction axis '$axis'. Valid axes are :consumer, :prey, or an existing population subgroup symbol.",
            ),
        )
    end
