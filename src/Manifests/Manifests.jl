"""Model setup export, import, and replay."""
module Manifests

using Dates
using JSON

include("serialization.jl")

export export_manifest, construct_from_manifest

const MODEL_SETUP_SCHEMA = "agate.model_setup.v1"
const MODEL_SETUP_KEYS = ("schema", "created_at", "agate", "model", "kwargs")
const MODEL_KEYS = ("family",)

function check_keys(x, allowed, path)
    x isa AbstractDict || throw(ArgumentError("$path must be an object."))
    for key in keys(x)
        key in allowed || throw(ArgumentError("$path has unsupported field $(repr(key))."))
    end
    return x
end

function required(x::AbstractDict, key, path)
    haskey(x, key) || throw(ArgumentError("$path is missing required field $(repr(key))."))
    return x[key]
end

agate_version() = string(something(Base.pkgversion(parentmodule(@__MODULE__)), "unknown"))

"""
    export_manifest(path, setup) -> path

Write an Agate model setup manifest to `path` as pretty-printed JSON.
"""
function export_manifest(path::AbstractString, setup::AbstractDict)
    open(path, "w") do io
        JSON.print(io, setup, 4)
        println(io)
    end
    return path
end

function export_manifest(::AbstractString, setup)
    throw(ArgumentError("Expected an Agate model setup dictionary, got $(typeof(setup))."))
end

"""
    construct_from_manifest(setup; grid=nothing, arch=nothing) -> bgc
    construct_from_manifest(path::AbstractString; grid=nothing, arch=nothing) -> bgc

Reconstruct an Agate biogeochemistry object from an exported model setup.
"""
function construct_from_manifest(setup::AbstractDict; grid=nothing, arch=nothing)
    family = model_family(setup)
    kwargs = constructor_kwargs(setup, family; grid, arch)
    models = getfield(parentmodule(@__MODULE__), :Models)
    return getfield(models, Symbol(family)).construct(; kwargs...)
end

construct_from_manifest(path::AbstractString; grid=nothing, arch=nothing) =
    construct_from_manifest(JSON.parsefile(path); grid, arch)

function manifest_group_entries(groups, diameters_by_group)
    return Any[
        Dict{String,Any}(
            "name" => string(group),
            "diameters" => diameters_by_group[string(group)],
        ) for group in groups
    ]
end

function manifest_size_structure(group_roles, diameters_by_group)
    return Dict{String,Any}(
        "phytoplankton" => manifest_group_entries(
            group_roles.phytoplankton, diameters_by_group
        ),
        "zooplankton" => manifest_group_entries(
            group_roles.zooplankton, diameters_by_group
        ),
    )
end

function default_model_manifest(family::Symbol, data; group_roles=nothing)
    family_name = string(family)
    kwargs = Dict{String,Any}(
        "parameters" => data.parameter_values,
        "sinking_tracers" => data.sinking_tracers,
        "open_bottom" => data.open_bottom,
        "scalar_type" => data.scalar_type,
    )

    if isnothing(group_roles)
        kwargs["phyto_size_structure"] = data.plankton_diameters_by_group["P"]
        kwargs["zoo_size_structure"] = data.plankton_diameters_by_group["Z"]
    else
        kwargs["size_structure"] =
            manifest_size_structure(group_roles, data.plankton_diameters_by_group)
    end

    return Dict{String,Any}(
        "schema" => MODEL_SETUP_SCHEMA,
        "created_at" => string(now(UTC)),
        "agate" => Dict{String,Any}(
            "version" => agate_version(),
            "julia_version" => string(VERSION),
        ),
        "model" => Dict{String,Any}("family" => family_name),
        "kwargs" => kwargs,
    )
end

function model_family(setup::AbstractDict)
    check_keys(setup, MODEL_SETUP_KEYS, "Model setup")
    schema = required(setup, "schema", "Model setup")
    schema == MODEL_SETUP_SCHEMA || throw(
        ArgumentError(
            "Unsupported Agate model setup schema $(repr(schema)); expected $(repr(MODEL_SETUP_SCHEMA))."
        ),
    )

    model = check_keys(required(setup, "model", "Model setup"), MODEL_KEYS, "Model setup model")
    family = required(model, "family", "Model setup model")
    family in ("DARWIN", "NiPiZD") && return family

    throw(ArgumentError("Unsupported model family $(repr(family))."))
end

function constructor_kwargs(setup::AbstractDict, family; grid=nothing, arch=nothing)
    common = ("parameters", "sinking_tracers", "open_bottom", "scalar_type")
    allowed = family == "NiPiZD" ? (common..., "size_structure") :
        (common..., "phyto_size_structure", "zoo_size_structure")
    kwargs = check_keys(required(setup, "kwargs", "Model setup"), allowed, "Model setup kwargs")
    for key in allowed
        required(kwargs, key, "Model setup kwargs")
    end

    pairs = Pair{Symbol,Any}[]
    if family == "NiPiZD"
        push!(pairs, :size_structure => named_size_structure_kwargs(kwargs["size_structure"]))
    else
        for key in ("phyto_size_structure", "zoo_size_structure")
            value = kwargs[key]
            value isa AbstractVector ||
                throw(ArgumentError("Model setup kwargs.$key must be an array."))
            push!(pairs, Symbol(key) => setup_value(value))
        end
    end

    open_bottom = kwargs["open_bottom"]
    open_bottom isa Bool ||
        throw(ArgumentError("Model setup kwargs.open_bottom must be a boolean."))
    push!(pairs, :open_bottom => open_bottom)
    push!(pairs, :parameters => parameter_kwargs(kwargs["parameters"]))
    push!(pairs, :sinking_tracers => sinking_tracers_kwargs(kwargs["sinking_tracers"]))
    push!(pairs, :scalar_type => decode_scalar_type(kwargs["scalar_type"]))

    !isnothing(grid) && push!(pairs, :grid => grid)
    !isnothing(arch) && push!(pairs, :arch => arch)

    return (; pairs...)
end

function parameter_kwargs(parameters)
    parameters isa AbstractDict ||
        throw(ArgumentError("Model setup kwargs.parameters must be an object."))
    return (; (Symbol(k) => parameter_value(v, k) for (k, v) in pairs(parameters))...)
end

function named_size_structure_kwargs(size_structure)
    path = "Model setup kwargs.size_structure"
    check_keys(size_structure, ("phytoplankton", "zooplankton"), path)
    return (;
        phytoplankton=named_role_size_structure(size_structure, "phytoplankton"),
        zooplankton=named_role_size_structure(size_structure, "zooplankton"),
    )
end

function named_role_size_structure(size_structure::AbstractDict, role::AbstractString)
    path = "Model setup kwargs.size_structure.$role"
    groups = required(size_structure, role, "Model setup kwargs.size_structure")
    groups isa AbstractVector || throw(ArgumentError("$path must be an array."))

    entries = Pair{Symbol,Any}[]
    for (i, group) in pairs(groups)
        group_path = "$path[$i]"
        group = check_keys(group, ("name", "diameters"), group_path)
        name = required(group, "name", group_path)
        name isa AbstractString || throw(ArgumentError("$group_path.name must be a string."))
        diameters = required(group, "diameters", group_path)
        diameters isa AbstractVector ||
            throw(ArgumentError("$group_path.diameters must be an array."))
        push!(entries, Symbol(name) => setup_value(diameters))
    end
    return (; entries...)
end

function sinking_tracers_kwargs(sinking)
    isnothing(sinking) && return nothing
    path = "Model setup kwargs.sinking_tracers"
    sinking isa AbstractVector || throw(ArgumentError("$path must be an array or null."))

    entries = Pair{Symbol,Any}[]
    for (i, item) in pairs(sinking)
        item_path = "$path[$i]"
        item = check_keys(item, ("name", "value"), item_path)
        name = required(item, "name", item_path)
        name isa AbstractString || throw(ArgumentError("$item_path.name must be a string."))
        push!(entries, Symbol(name) => setup_value(required(item, "value", item_path)))
    end
    return (; entries...)
end

function parameter_value(value, name)
    if value isa AbstractVector && all(row -> row isa AbstractVector, value)
        rows = Any[row for row in value]
        isempty(rows) && return Matrix{Any}(undef, 0, 0)
        ncols = length(rows[1])
        all(row -> length(row) == ncols, rows) || throw(
            ArgumentError("Model setup parameter $(repr(name)) matrix must be rectangular."),
        )
        return [setup_value(rows[i][j]) for i in eachindex(rows), j in eachindex(rows[1])]
    end
    return setup_value(value)
end

setup_value(x) = x
function setup_value(x::AbstractString)
    x == "NaN" && return NaN
    x == "Inf" && return Inf
    x == "-Inf" && return -Inf
    return x
end
setup_value(x::AbstractVector) = Any[setup_value(v) for v in x]
setup_value(x::AbstractDict) = (; (Symbol(k) => setup_value(v) for (k, v) in pairs(x))...)

function decode_scalar_type(x)
    isnothing(x) && return nothing
    x isa Type && x <: Real && return x
    x isa AbstractString || throw(ArgumentError("Model setup scalar_type must be a string, got $(typeof(x))."))

    name = Symbol(x)
    for mod in (Core, Base)
        if isdefined(mod, name)
            T = getfield(mod, name)
            T isa Type && T <: Real && return T
        end
    end

    throw(ArgumentError("Unsupported model setup scalar_type $(repr(x))."))
end

end
