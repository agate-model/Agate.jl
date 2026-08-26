"""Return parameter specs stored on the canonical consumer-by-prey axes."""
function interaction_parameter_specs(source)
    return Tuple(
        name => spec for (name, spec) in pairs(parameter_directory(source)) if
        spec.axes == (:consumer, :prey)
    )
end

"""Return runtime metadata for consumer-by-prey parameter matrices."""
function interaction_axis_metadata(source, layout::ModelLayout)
    specs = interaction_parameter_specs(source)
    isempty(specs) && return nothing

    consumer_classes = Tuple(layout.class_symbols[i] for i in layout.consumer_indices)
    prey_classes = Tuple(layout.class_symbols[i] for i in layout.prey_indices)
    parameter_names = Tuple(first(entry) for entry in specs)

    return (;
        parameters=parameter_names,
        consumers=consumer_classes,
        prey=prey_classes,
    )
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
