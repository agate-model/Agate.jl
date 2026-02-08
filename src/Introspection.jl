"""Introspection helpers for constructed Agate biogeochemistry instances.

These functions are meant for interactive use (REPL / notebooks): they let you
inspect what a constructed model exposes without digging into generated types.

All helpers are CPU/GPU-safe: they operate on small metadata (names and keys)
and never touch state arrays.
"""
module Introspection

export tracer_names
export auxiliary_field_names
export parameter_names
export model_summary
export describe

import Oceananigans.Biogeochemistry:
    required_biogeochemical_auxiliary_fields, required_biogeochemical_tracers

@inline function _preview_list(xs; n::Int=12)
    m = length(xs)
    if m <= n
        return join(string.(xs), ", ")
    end
    head = join(string.(xs[1:n]), ", ")
    return head * ", ... (" * string(m) * ")"
end

"""
    tracer_names(bgc) -> Vector{Symbol}

Return the ordered tracer symbols required by `bgc`.

The ordering matches Oceananigans / OceanBioME state-vector conventions.
"""
function tracer_names(bgc)::Vector{Symbol}
    return collect(required_biogeochemical_tracers(bgc))
end

"""
    auxiliary_field_names(bgc) -> Vector{Symbol}

Return the ordered auxiliary field symbols required by `bgc`.

Auxiliary fields are non-tracer state fields (for example, light or temperature)
that appear in tracer tendencies.
"""
function auxiliary_field_names(bgc)::Vector{Symbol}
    return collect(required_biogeochemical_auxiliary_fields(bgc))
end

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

"""
    model_summary(bgc) -> NamedTuple

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
        verbose ? join(string.(s.tracers), ", ") : _preview_list(s.tracers),
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
                _preview_list(s.auxiliary_fields)
            end,
        )
    end

    println(
        io,
        "  parameters (",
        length(s.parameters),
        "): ",
        verbose ? join(string.(s.parameters), ", ") : _preview_list(s.parameters),
    )
    println(io, "  sinking velocities: ", s.has_sinking_velocities ? "yes" : "no")

    return nothing
end

describe(bgc; kwargs...) = describe(stdout, bgc; kwargs...)

end # module
