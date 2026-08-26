"""Introspection helpers for constructed Agate biogeochemistry instances.

These functions are meant for interactive use (REPL / notebooks): they let you
inspect what a constructed model exposes without digging into generated types.
"""
module Introspection

export tracer_names
export auxiliary_field_names
export parameter_names
export plankton_groups
export plankton_tracers
export plankton_diameters
export nonplankton_tracers
export tracer_groups
export interaction_matrix
export model_summary
export describe

import Oceananigans.Biogeochemistry:
    required_biogeochemical_auxiliary_fields, required_biogeochemical_tracers


@inline function preview_list(xs; n::Int=12)
    m = length(xs)
    if m <= n
        return join(string.(xs), ", ")
    end
    head = join(string.(xs[1:n]), ", ")
    return head * ", ... (" * string(m) * ")"
end

"""    tracer_names(bgc) -> Vector{Symbol}

Return the ordered tracer symbols required by `bgc`.

This helper is intended for interactive inspection, so it materializes the
underlying tracer-name tuple as a `Vector{Symbol}`.

The ordering matches Oceananigans / OceanBioME state-vector conventions.
"""
@inline tracer_names(bgc)::Vector{Symbol} = collect(required_biogeochemical_tracers(bgc))

"""    auxiliary_field_names(bgc) -> Vector{Symbol}

Return the ordered auxiliary field symbols required by `bgc`.

Auxiliary fields are non-tracer state fields (for example, light or temperature)
that appear in tracer tendencies.
"""
@inline auxiliary_field_names(bgc)::Vector{Symbol} = collect(
    required_biogeochemical_auxiliary_fields(bgc)
)

"""
    parameter_names(bgc) -> Vector{Symbol}

Return the parameter keys available on `bgc.parameters`.

This list describes the resolved parameter fields available on the constructed
biogeochemistry instance.
"""
function parameter_names(bgc)::Vector{Symbol}
    params = getproperty(bgc, :parameters)
    return collect(propertynames(params))
end

function _model_metadata(bgc)
    hasproperty(bgc, :metadata) || return nothing
    return getproperty(bgc, :metadata)
end

"""    plankton_groups(bgc) -> NamedTuple

Return a `NamedTuple` mapping plankton group symbols to ecological class symbols.
For multi-state populations each size class appears once, independent of the number
of physical prognostic state tracers. Group order follows the realized model layout.
"""
function plankton_groups(bgc)
    metadata = _model_metadata(bgc)
    metadata === nothing && return NamedTuple()
    names = keys(metadata.group_classes)
    return NamedTuple{names}(
        Tuple(collect(classes) for classes in values(metadata.group_classes))
    )
end

"""    plankton_tracers(bgc) -> Vector{Symbol}

Return all plankton tracer symbols as a flat vector in runtime group order.
"""
function plankton_tracers(bgc)
    metadata = _model_metadata(bgc)
    metadata === nothing && return Symbol[]
    return collect(metadata.population_tracers)
end

"""    plankton_diameters(bgc) -> Vector

Return the equivalent spherical diameters for realized plankton ecological classes.
The ordering follows the flattened values of `plankton_groups(bgc)`, with one diameter
per size class even when each class carries multiple prognostic state tracers. Models
without plankton diameter metadata return an empty vector.
"""
function plankton_diameters(bgc)
    metadata = _model_metadata(bgc)
    metadata === nothing && return []
    return collect(metadata.plankton_diameters)
end

"""    nonplankton_tracers(bgc) -> Vector{Symbol}

Return the tracer symbols that are not part of a plankton group.
"""
function nonplankton_tracers(bgc)
    plankton = Set(plankton_tracers(bgc))
    return [tracer for tracer in tracer_names(bgc) if tracer ∉ plankton]
end

"""    tracer_groups(bgc) -> NamedTuple

Return a structural grouping summary of the constructed tracer layout.
"""
function tracer_groups(bgc)
    return (
        all=tracer_names(bgc),
        plankton=plankton_tracers(bgc),
        nonplankton=nonplankton_tracers(bgc),
        by_group=plankton_groups(bgc),
    )
end

function _interaction_axes(bgc)
    metadata = _model_metadata(bgc)
    metadata === nothing && throw(ArgumentError("No interaction matrix metadata found for this model."))
    axes = metadata.interaction_axes
    axes === nothing && throw(ArgumentError("No interaction matrices found for this model."))
    return axes
end

function _require_interaction_parameter(bgc, kind::Symbol, axes)
    kind in axes.parameters || begin
        available_text = isempty(axes.parameters) ? "none" : join(string.(axes.parameters), ", ")
        throw(
            ArgumentError(
                "Unknown interaction matrix parameter: $kind. Available parameters are: $available_text."
            ),
        )
    end
    hasproperty(bgc.parameters, kind) || throw(
        ArgumentError("Interaction parameter :$kind is missing from runtime parameters."),
    )
    matrix = getproperty(bgc.parameters, kind)
    matrix isa AbstractMatrix || throw(
        ArgumentError("Interaction parameter :$kind is not stored as a matrix."),
    )
    return matrix
end

function _require_interaction_shape(matrix, rows, columns, kind::Symbol)
    expected = (length(rows), length(columns))
    size(matrix) == expected && return nothing

    throw(
        ArgumentError(
            "Interaction matrix $kind has size $(size(matrix)); expected $expected from labelled interaction axes.",
        ),
    )
end

"""    interaction_matrix(bgc, parameter::Symbol) -> NamedTuple

Return a consumer-by-prey parameter matrix with ecological class labels.

`parameter` is the canonical model parameter identity, for example
`:palatability_matrix` or `:assimilation_matrix`. The returned `NamedTuple` contains
`parameter`, `matrix`, `rows`, `columns`, `row_axis`, and `column_axis`.
"""
function interaction_matrix(bgc, parameter::Symbol)
    axes = _interaction_axes(bgc)
    matrix = _require_interaction_parameter(bgc, parameter, axes)
    rows = collect(axes.consumers)
    columns = collect(axes.prey)
    _require_interaction_shape(matrix, rows, columns, parameter)

    return (
        parameter=parameter,
        matrix=matrix,
        rows=rows,
        columns=columns,
        row_axis=:consumer,
        column_axis=:prey,
    )
end

"""    model_summary(bgc) -> NamedTuple

Return a compact summary of a constructed biogeochemistry instance.

The returned `NamedTuple` contains:

- `tracers::Vector{Symbol}`
- `auxiliary_fields::Vector{Symbol}`
- `parameters::Vector{Symbol}`
- `has_sinking_velocities::Bool`
"""
function model_summary(bgc)
    return (
        tracers=tracer_names(bgc),
        auxiliary_fields=auxiliary_field_names(bgc),
        parameters=parameter_names(bgc),
        has_sinking_velocities=Base.hasproperty(bgc, :sinking_velocities) && getproperty(bgc, :sinking_velocities) !== nothing,
    )
end

"""
    describe([io], bgc; verbose=false)

Print a human-readable summary of `bgc`.

Set `verbose=true` to print full tracer / parameter lists.
"""
function describe(io::IO, bgc; verbose::Bool=false)
    s = model_summary(bgc)

    println(io, "Agate biogeochemistry instance")
    println(
        io,
        "  tracers (",
        length(s.tracers),
        "): ",
        verbose ? join(string.(s.tracers), ", ") : preview_list(s.tracers),
    )

    if isempty(s.auxiliary_fields)
        println(io, "  auxiliary fields: (none)")
    else
        println(
            io,
            "  auxiliary fields (",
            length(s.auxiliary_fields),
            "): ",
            if verbose
                join(string.(s.auxiliary_fields), ", ")
            else
                preview_list(s.auxiliary_fields)
            end,
        )
    end

    println(
        io,
        "  parameters (",
        length(s.parameters),
        "): ",
        verbose ? join(string.(s.parameters), ", ") : preview_list(s.parameters),
    )
    println(io, "  sinking velocities: ", s.has_sinking_velocities ? "yes" : "no")

    return nothing
end

describe(bgc; kwargs...) = describe(stdout, bgc; kwargs...)

end # module
