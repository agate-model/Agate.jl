"""Return runtime metadata for parameters stored on canonical consumer-by-prey axes."""
function interaction_axis_metadata(plan, layout::ModelLayout)
    consumers = Tuple(layout.class_symbols[index] for index in layout.consumer_indices)
    prey = Tuple(layout.class_symbols[index] for index in layout.prey_indices)
    names = Tuple(
        name for name in keys(plan.parameters)
        if begin
            parameter = getproperty(plan.parameters, name)
            parameter.runtime_bound &&
                parameter.axes == (:consumer, :resource) &&
                parameter.storage_labels == (consumers, prey)
        end
    )
    isempty(names) && return nothing
    return (; parameters=names, consumers, prey)
end
