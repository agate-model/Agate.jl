"""Model setup export, import, and replay."""
module Manifests

using Dates
using JSON
using TOML

include("serialization.jl")

export export_manifest, construct_from_manifest

const MODEL_SETUP_SCHEMA = "agate.model_setup.v1"

agate_project_root() = dirname(dirname(@__DIR__))
agate_project_toml() = joinpath(agate_project_root(), "Project.toml")

function agate_version()
    project = TOML.parsefile(agate_project_toml())
    return string(project["version"])
end

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
    kwargs = constructor_kwargs(setup; grid, arch)
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
    get(setup, "schema", nothing) == MODEL_SETUP_SCHEMA ||
        throw(ArgumentError("Unsupported Agate model setup schema."))

    model = get(setup, "model", nothing)
    model isa AbstractDict || throw(ArgumentError("Model setup is missing a model section."))

    family = get(model, "family", nothing)
    family in ("DARWIN", "NiPiZD") && return family

    throw(ArgumentError("Unsupported model family $(repr(family))."))
end

function constructor_kwargs(setup::AbstractDict; grid=nothing, arch=nothing)
    kwargs = get(setup, "kwargs", nothing)
    kwargs isa AbstractDict || throw(ArgumentError("Model setup is missing constructor kwargs."))

    pairs = Pair{Symbol,Any}[]

    haskey(kwargs, "size_structure") &&
        push!(pairs, :size_structure => named_size_structure_kwargs(kwargs["size_structure"]))
    for key in ("phyto_size_structure", "zoo_size_structure", "open_bottom")
        haskey(kwargs, key) && push!(pairs, Symbol(key) => setup_value(kwargs[key]))
    end

    haskey(kwargs, "parameters") &&
        push!(pairs, :parameters => parameter_kwargs(kwargs["parameters"]))
    haskey(kwargs, "sinking_tracers") &&
        push!(pairs, :sinking_tracers => sinking_tracers_kwargs(kwargs["sinking_tracers"]))
    haskey(kwargs, "scalar_type") &&
        push!(pairs, :scalar_type => decode_scalar_type(kwargs["scalar_type"]))

    !isnothing(grid) && push!(pairs, :grid => grid)
    !isnothing(arch) && push!(pairs, :arch => arch)

    return (; pairs...)
end

function parameter_kwargs(parameters::AbstractDict)
    return (; (Symbol(k) => parameter_value(v) for (k, v) in pairs(parameters))...)
end

function named_size_structure_kwargs(size_structure::AbstractDict)
    return (;
        phytoplankton=named_role_size_structure(size_structure, "phytoplankton"),
        zooplankton=named_role_size_structure(size_structure, "zooplankton"),
    )
end

function named_size_structure_kwargs(size_structure)
    throw(ArgumentError("Model setup size_structure must be an object."))
end

function named_role_size_structure(size_structure::AbstractDict, role::AbstractString)
    groups = get(size_structure, role, nothing)
    groups isa AbstractVector ||
        throw(ArgumentError("Model setup size_structure.$role must be an array."))
    return (;
        (Symbol(group["name"]) => setup_value(group["diameters"]) for group in groups)...
    )
end

function sinking_tracers_kwargs(sinking)
    isnothing(sinking) && return nothing
    return (; (Symbol(item["name"]) => setup_value(item["value"]) for item in sinking)...)
end

function parameter_value(value)
    if value isa AbstractVector && all(row -> row isa AbstractVector, value)
        rows = Any[row for row in value]
        isempty(rows) && return Matrix{Any}(undef, 0, 0)
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
