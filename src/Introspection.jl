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

using ..Runtime: TracerIndex

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
    keys = collect(propertynames(params))

    # Some parameter fields are internal containers (for example, interaction
    # matrices plus axis maps). Hide these from the public parameter list.
    filter!(k -> k !== :interactions, keys)
    return keys
end

function _tracer_index(bgc)
    hasproperty(bgc, :tracers) || return nothing
    tracers = getproperty(bgc, :tracers)
    hasproperty(tracers, :idx) || return nothing
    return getproperty(tracers, :idx)
end

_all_tracer_symbols(::TracerIndex{TR,GS,AF,NG}) where {TR,GS,AF,NG} = collect(TR)
_plankton_group_symbols(::TracerIndex{TR,GS,AF,NG}) where {TR,GS,AF,NG} = collect(GS)

"""    plankton_groups(bgc) -> NamedTuple

Return a `NamedTuple` mapping plankton group symbols to the tracer symbols in
each group. The group and tracer order follows the constructed runtime tracer
layout.
"""
function plankton_groups(bgc)
    idx = _tracer_index(bgc)
    idx isa TracerIndex || return NamedTuple()

    all_tracers = _all_tracer_symbols(idx)
    groups = _plankton_group_symbols(idx)
    isempty(groups) && return NamedTuple()

    pairs = Pair{Symbol,Vector{Symbol}}[]
    for (i, group) in enumerate(groups)
        first_index = idx.group_bases[i]
        n_tracers = idx.group_counts[i]
        push!(pairs, group => all_tracers[first_index:(first_index + n_tracers - 1)])
    end

    return (; pairs...)
end

"""    plankton_tracers(bgc) -> Vector{Symbol}

Return all plankton tracer symbols as a flat vector in runtime group order.
"""
function plankton_tracers(bgc)
    groups = plankton_groups(bgc)
    isempty(groups) && return Symbol[]

    tracers = Symbol[]
    for group_tracers in values(groups)
        append!(tracers, group_tracers)
    end

    return tracers
end

"""    plankton_diameters(bgc) -> Vector

Return the equivalent spherical diameters for plankton tracers in the same order
as `plankton_tracers(bgc)`. Models without plankton diameter metadata return
an empty vector.
"""
function plankton_diameters(bgc)
    hasproperty(bgc, :plankton_diameters) || return []
    return collect(getproperty(bgc, :plankton_diameters))
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

function _interaction_container(bgc)
    hasproperty(bgc, :parameters) || return nothing
    params = getproperty(bgc, :parameters)
    hasproperty(params, :interactions) || return nothing
    return getproperty(params, :interactions)
end

const _INTERACTION_AXIS_FIELDS = (
    :consumer_global,
    :prey_global,
    :global_to_consumer,
    :global_to_prey,
)


function _available_interaction_kinds(interactions)
    kinds = Symbol[]

    for property in propertynames(interactions)
        property in _INTERACTION_AXIS_FIELDS && continue
        value = getproperty(interactions, property)
        value isa AbstractMatrix && push!(kinds, property)
    end

    return kinds
end

function _require_interactions(bgc)
    interactions = _interaction_container(bgc)
    interactions === nothing &&
        throw(ArgumentError("No interaction matrices found for this model."))
    return interactions
end

function _require_interaction_kind(interactions, kind::Symbol)
    available = _available_interaction_kinds(interactions)
    kind in available && return getproperty(interactions, kind)

    available_text = isempty(available) ? "none" : join(string.(available), ", ")
    throw(
        ArgumentError(
            "Unknown interaction matrix kind: $kind. Available kinds are: $available_text."
        ),
    )
end

function _plankton_axis_labels(bgc, indices)
    plankton = plankton_tracers(bgc)
    isempty(plankton) &&
        throw(ArgumentError("Interaction axes require plankton tracer metadata."))

    labels = Symbol[]
    for index in Array(indices)
        i = Int(index)
        1 <= i <= length(plankton) || throw(
            ArgumentError(
                "Interaction axis index $i is outside the plankton tracer axis of length $(length(plankton)).",
            ),
        )
        push!(labels, plankton[i])
    end

    return labels
end

function _interaction_axes(bgc, interactions)
    for field in (:consumer_global, :prey_global)
        hasproperty(interactions, field) || throw(
            ArgumentError("Interaction matrices are missing required axis field: $field."),
        )
    end

    rows = _plankton_axis_labels(bgc, getproperty(interactions, :consumer_global))
    columns = _plankton_axis_labels(bgc, getproperty(interactions, :prey_global))

    return (rows=rows, columns=columns, row_axis=:consumer, column_axis=:prey)
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

"""    interaction_matrix(bgc, kind::Symbol) -> NamedTuple

Return an interaction matrix with consumer and prey labels.

Supported `kind` values are the available interaction matrix fields, such as
`:palatability` and `:assimilation`. The returned `NamedTuple` contains `kind`,
`matrix`, `rows`, `columns`, `row_axis`, and `column_axis`. Matrix orientation
follows the runtime container: rows are consumers and columns are prey.
"""
function interaction_matrix(bgc, kind::Symbol)
    interactions = _require_interactions(bgc)
    axes = _interaction_axes(bgc, interactions)
    matrix = _require_interaction_kind(interactions, kind)
    _require_interaction_shape(matrix, axes.rows, axes.columns, kind)

    return (
        kind=kind,
        matrix=matrix,
        rows=axes.rows,
        columns=axes.columns,
        row_axis=axes.row_axis,
        column_axis=axes.column_axis,
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
        has_sinking_velocities=Base.hasproperty(bgc, :sinking_velocities),
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
