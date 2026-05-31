"""Construction manifest I/O, serialization, and replay."""
module Manifests

using Dates
using JSON
using TOML

include("serialization.jl")

export export_manifest, construct_from_manifest

const CONSTRUCTION_MANIFEST_SCHEMA = "agate.construction_manifest.v1"
const MODEL_CONSTRUCTOR_RECIPE_TYPE = "model_constructor"
const MODEL_CONSTRUCTOR_FAMILIES = ("DARWIN", "NiPiZD")

agate_project_root() = dirname(dirname(@__DIR__))
agate_project_toml() = joinpath(agate_project_root(), "Project.toml")

function agate_version()
    project = TOML.parsefile(agate_project_toml())
    return string(project["version"])
end

"""
    export_manifest(path, manifest) -> path

Write an Agate construction manifest dictionary to `path` as pretty-printed JSON.
"""
function export_manifest(path::AbstractString, manifest::AbstractDict)
    open(path, "w") do io
        JSON.print(io, manifest, 4)
        println(io)
    end
    return path
end

function export_manifest(::AbstractString, manifest)
    throw(ArgumentError("Expected an Agate construction manifest dictionary, got $(typeof(manifest))."))
end

"""
    construct_from_manifest(manifest; grid=nothing, arch=nothing) -> bgc
    construct_from_manifest(path::AbstractString; grid=nothing, arch=nothing) -> bgc

Reconstruct an Agate biogeochemistry object from a construction manifest.
"""
function construct_from_manifest(manifest::AbstractDict; grid=nothing, arch=nothing)
    recipe = validated_recipe(manifest)
    family = model_constructor_family(recipe)
    kwargs = recipe_constructor_kwargs(recipe["kwargs"], manifest; grid=grid, arch=arch)
    models = getfield(parentmodule(@__MODULE__), :Models)
    return getfield(models, Symbol(family)).construct(; kwargs...)
end

construct_from_manifest(path::AbstractString; grid=nothing, arch=nothing) =
    construct_from_manifest(JSON.parsefile(path); grid=grid, arch=arch)

function default_model_recipe(family::Symbol, data)
    family_name = string(family)
    diameters = data.plankton_diameters_by_group
    return Dict{String,Any}(
        "type" => MODEL_CONSTRUCTOR_RECIPE_TYPE,
        "family" => family_name,
        "kwargs" => Dict{String,Any}(
            "phyto_size_structure" => diameters["P"],
            "zoo_size_structure" => diameters["Z"],
            "parameters" => data.parameter_values,
            "sinking_tracers" => data.sinking_tracers_recipe,
            "open_bottom" => data.sinking["open_bottom"],
            "scalar_type" => data.scalar_type,
        ),
    )
end

function default_model_manifest(family::Symbol, data)
    family_name = string(family)
    recipe = default_model_recipe(family, data)
    return Dict{String,Any}(
        "schema" => CONSTRUCTION_MANIFEST_SCHEMA,
        "created_at" => string(now(UTC)),
        "agate" => Dict{String,Any}(
            "version" => agate_version(),
            "julia_version" => string(VERSION),
        ),
        "model" => Dict{String,Any}(
            "id" => string(family_name, "/default"),
            "family" => family_name,
            "citation" => "default",
            "tag" => "default",
        ),
        "recipe" => recipe,
        "resolved" => Dict{String,Any}(
            "tracers" => data.tracers,
            "auxiliary_fields" => data.auxiliary_fields,
            "parameters" => data.parameters,
            "plankton_diameters" => data.plankton_diameters,
            "plankton_diameters_by_group" => data.plankton_diameters_by_group,
            "scalar_type" => data.scalar_type,
            "architecture" => data.architecture,
            "has_sinking_velocities" => data.has_sinking_velocities,
            "sinking" => data.sinking,
        ),
    )
end

function validated_recipe(manifest::AbstractDict)
    schema = get(manifest, "schema", nothing)
    schema isa AbstractString || throw(ArgumentError("Manifest is missing a string \"schema\" entry."))
    schema == CONSTRUCTION_MANIFEST_SCHEMA || throw(
        ArgumentError("Unsupported manifest schema $(repr(schema)); expected $(repr(CONSTRUCTION_MANIFEST_SCHEMA)).")
    )

    recipe = get(manifest, "recipe", nothing)
    recipe isa AbstractDict || throw(
        ArgumentError("Manifest does not contain a dictionary \"recipe\" section required for reconstruction.")
    )

    recipe_type = get(recipe, "type", nothing)
    recipe_type isa AbstractString || throw(ArgumentError("Manifest recipe is missing a string \"type\" entry."))
    recipe_type == MODEL_CONSTRUCTOR_RECIPE_TYPE || throw(
        ArgumentError("Unsupported manifest recipe type $(repr(recipe_type)); expected $(repr(MODEL_CONSTRUCTOR_RECIPE_TYPE)).")
    )

    family = get(recipe, "family", nothing)
    family isa AbstractString || throw(ArgumentError("Manifest recipe is missing a string \"family\" entry."))
    model_constructor_family(recipe)

    kwargs = get(recipe, "kwargs", nothing)
    kwargs isa AbstractDict || throw(ArgumentError("Manifest recipe is missing a dictionary \"kwargs\" entry."))

    return recipe
end

function model_constructor_family(recipe::AbstractDict)
    family = recipe["family"]
    family in MODEL_CONSTRUCTOR_FAMILIES && return family
    throw(ArgumentError("Unsupported model constructor family $(repr(family))."))
end

function recipe_constructor_kwargs(kwargs_dict::AbstractDict, manifest::AbstractDict; grid=nothing, arch=nothing)
    pairs = Pair{Symbol,Any}[]

    for key in ("phyto_size_structure", "zoo_size_structure", "open_bottom")
        haskey(kwargs_dict, key) && push!(pairs, Symbol(key) => kwargs_dict[key])
    end

    if haskey(kwargs_dict, "parameters")
        parameters = kwargs_dict["parameters"]
        parameters isa AbstractDict || throw(
            ArgumentError("Manifest recipe \"parameters\" entry must be a dictionary, got $(typeof(parameters)).")
        )
        push!(pairs, :parameters => parameter_kwargs(parameters, manifest))
    end

    haskey(kwargs_dict, "sinking_tracers") &&
        push!(pairs, :sinking_tracers => sinking_tracers_kwargs(kwargs_dict["sinking_tracers"]))

    haskey(kwargs_dict, "scalar_type") &&
        push!(pairs, :scalar_type => decode_scalar_type(kwargs_dict["scalar_type"]))

    !isnothing(grid) && push!(pairs, :grid => grid)
    !isnothing(arch) && push!(pairs, :arch => arch)

    return (; pairs...)
end

function sinking_tracers_kwargs(sinking)
    isnothing(sinking) && return nothing
    if sinking isa AbstractDict
        return dict_to_namedtuple(sinking)
    elseif sinking isa AbstractVector
        pairs = Pair{Symbol,Any}[]
        for item in sinking
            item isa AbstractDict || throw(
                ArgumentError("Manifest recipe ordered sinking tracer entries must be dictionaries, got $(typeof(item)).")
            )
            haskey(item, "name") && haskey(item, "value") || throw(
                ArgumentError("Manifest recipe ordered sinking tracer entries must contain \"name\" and \"value\".")
            )
            push!(pairs, Symbol(item["name"]) => manifest_value_to_julia(item["value"]))
        end
        return (; pairs...)
    else
        throw(ArgumentError("Manifest recipe \"sinking_tracers\" entry must be an ordered vector, dictionary, or nothing, got $(typeof(sinking))."))
    end
end

function parameter_kwargs(parameters::AbstractDict, manifest::AbstractDict)
    records = get(get(manifest, "resolved", Dict{String,Any}()), "parameters", Dict{String,Any}())
    return (; (
        Symbol(k) => manifest_parameter_value_to_julia(k, v, records) for
        (k, v) in pairs(parameters)
    )...)
end

function manifest_parameter_value_to_julia(key, value, records)
    record = records isa AbstractDict ? get(records, string(key), nothing) : nothing
    if record isa AbstractDict && get(record, "shape", nothing) == "matrix"
        return manifest_matrix(value)
    end
    return manifest_value_to_julia(value)
end

function manifest_matrix(value)
    value isa AbstractVector || throw(
        ArgumentError("Manifest matrix parameter must be encoded as a vector of rows, got $(typeof(value)).")
    )
    rows = Any[value[i] for i in eachindex(value)]
    isempty(rows) && return Matrix{Any}(undef, 0, 0)
    all(row -> row isa AbstractVector, rows) || throw(ArgumentError("Manifest matrix parameter rows must be vectors."))
    ncols = length(rows[1])
    all(row -> length(row) == ncols, rows) || throw(
        ArgumentError("Manifest matrix parameter rows must all have the same length.")
    )
    return [manifest_value_to_julia(rows[i][j]) for i in eachindex(rows), j in 1:ncols]
end

dict_to_namedtuple(dict::AbstractDict) =
    (; (Symbol(k) => manifest_value_to_julia(v) for (k, v) in pairs(dict))...)

manifest_value_to_julia(x) = x

function decode_scalar_type(x)
    isnothing(x) && return nothing
    x isa Type && x <: Real && return x
    x isa AbstractString || throw(
        ArgumentError("Manifest recipe scalar_type must be a string or concrete Real type, got $(typeof(x)).")
    )

    name = Symbol(x)
    for mod in (Core, Base)
        if isdefined(mod, name)
            T = getfield(mod, name)
            T isa Type && T <: Real && return T
        end
    end

    throw(ArgumentError("Unsupported manifest scalar_type $(repr(x))."))
end

end
