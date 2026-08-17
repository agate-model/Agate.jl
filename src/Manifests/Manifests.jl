"""Model setup export, import, and replay."""
module Manifests

using Dates
using JSON

include("serialization.jl")

export export_manifest, construct_from_manifest, model_manifest

const MODEL_SETUP_SCHEMA = "agate.model_setup.v1"
# schema/model/kwargs and family kwargs are required for replay; created_at and agate
# are optional provenance when reading existing v1 files.
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
    return construct_from_manifest(Val(family), setup; grid, arch)
end

function construct_from_manifest(
    ::Val{family}, ::AbstractDict; grid=nothing, arch=nothing
) where {family}
    throw(ArgumentError("Unsupported model family $(repr(String(family)))."))
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

"""Build a resolved schema-v1 model manifest from shared data and family-specific kwargs."""
function model_manifest(family::Symbol, data, model_kwargs::AbstractDict)
    kwargs = Dict{String,Any}(
        "parameters" => data.parameter_values,
        "sinking_tracers" => data.sinking_tracers,
        "open_bottom" => data.open_bottom,
        "scalar_type" => data.scalar_type,
    )
    for (key, value) in pairs(model_kwargs)
        kwargs[string(key)] = value
    end

    return Dict{String,Any}(
        "schema" => MODEL_SETUP_SCHEMA,
        "created_at" => string(now(UTC)),
        "agate" => Dict{String,Any}(
            "version" => agate_version(),
            "julia_version" => string(VERSION),
        ),
        "model" => Dict{String,Any}("family" => string(family)),
        "kwargs" => kwargs,
    )
end

function default_model_manifest(family::Symbol, data; group_roles)
    model_kwargs = Dict{String,Any}(
        "size_structure" =>
            manifest_size_structure(group_roles, data.plankton_diameters_by_group),
    )
    return model_manifest(family, data, model_kwargs)
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
    family isa AbstractString ||
        throw(ArgumentError("Model setup model.family must be a string."))
    return Symbol(family)
end

const COMMON_CONSTRUCTOR_KEYS = (
    "parameters", "sinking_tracers", "open_bottom", "scalar_type"
)

function manifest_kwargs(setup::AbstractDict, model_keys::Tuple)
    allowed = (COMMON_CONSTRUCTOR_KEYS..., model_keys...)
    kwargs = check_keys(required(setup, "kwargs", "Model setup"), allowed, "Model setup kwargs")
    for key in allowed
        required(kwargs, key, "Model setup kwargs")
    end
    return kwargs
end

function common_constructor_kwargs(kwargs::AbstractDict; grid=nothing, arch=nothing)
    open_bottom = kwargs["open_bottom"]
    open_bottom isa Bool ||
        throw(ArgumentError("Model setup kwargs.open_bottom must be a boolean."))

    kwarg_pairs = Pair{Symbol,Any}[
        :open_bottom => open_bottom,
        :parameters => parameter_kwargs(kwargs["parameters"]),
        :sinking_tracers => sinking_tracers_kwargs(kwargs["sinking_tracers"]),
        :scalar_type => decode_scalar_type(kwargs["scalar_type"]),
    ]
    !isnothing(grid) && push!(kwarg_pairs, :grid => grid)
    !isnothing(arch) && push!(kwarg_pairs, :arch => arch)
    return (; kwarg_pairs...)
end

function size_structure_vector(kwargs::AbstractDict, key::AbstractString)
    value = kwargs[key]
    value isa AbstractVector ||
        throw(ArgumentError("Model setup kwargs.$key must be an array."))
    return setup_value(value)
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
    (value isa AbstractVector && all(row -> row isa AbstractVector, value)) ||
        return setup_value(value)
    isempty(value) && throw(
        ArgumentError(
            "Model setup parameter $(repr(name)) is an ambiguous empty array in schema v1."
        ),
    )
    ncols = length(first(value))
    all(row -> length(row) == ncols, value) || throw(
        ArgumentError("Model setup parameter $(repr(name)) matrix must be rectangular."),
    )
    return [setup_value(value[i][j]) for i in eachindex(value), j in 1:ncols]
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
