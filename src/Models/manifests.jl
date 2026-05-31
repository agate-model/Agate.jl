using Dates
using JSON
using TOML

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
    throw(
        ArgumentError(
            "Expected an Agate construction manifest dictionary, got $(typeof(manifest))."
        ),
    )
end



"""
    construct_from_manifest(manifest; grid=nothing, arch=nothing) -> bgc
    construct_from_manifest(path::AbstractString; grid=nothing, arch=nothing) -> bgc

Reconstruct an Agate biogeochemistry object from a model-level construction manifest.

This uses the manifest's `recipe` section to resolve an allowed model family.
The manifest must therefore contain a top-level `"recipe"` entry,
which is currently produced by model-level `construct_with_manifest` methods such as
`Agate.Models.DARWIN.construct_with_manifest` and
`Agate.Models.NiPiZD.construct_with_manifest`.

The manifest does not currently encode a grid specification. Pass `grid=...` when you
want to reconstruct against a specific grid; otherwise the target constructor's default
behavior is used.
"""
const CONSTRUCTION_MANIFEST_SCHEMA = "agate.construction_manifest.v1"
const MODEL_CONSTRUCTOR_RECIPE_TYPE = "model_constructor"
const MODEL_CONSTRUCTOR_FAMILIES = ("DARWIN", "NiPiZD")

function construct_from_manifest(manifest::AbstractDict; grid=nothing, arch=nothing)
    recipe = _validated_manifest_recipe(manifest)

    family = _model_constructor_family(recipe)
    kwargs_dict = recipe["kwargs"]
    constructor_kwargs = _recipe_constructor_kwargs(kwargs_dict, manifest; grid=grid, arch=arch)

    if family == "DARWIN"
        return DARWIN.construct(; constructor_kwargs...)
    elseif family == "NiPiZD"
        return NiPiZD.construct(; constructor_kwargs...)
    else
        throw(ArgumentError("Unsupported model constructor family $(repr(family))."))
    end
end

function _validated_manifest_recipe(manifest::AbstractDict)
    schema = get(manifest, "schema", nothing)
    schema isa AbstractString || throw(
        ArgumentError("Manifest is missing a string \"schema\" entry.")
    )
    schema == CONSTRUCTION_MANIFEST_SCHEMA || throw(
        ArgumentError(
            "Unsupported manifest schema $(repr(schema)); expected $(repr(CONSTRUCTION_MANIFEST_SCHEMA))."
        ),
    )

    if !haskey(manifest, "recipe")
        replayable = get(manifest, "replayable", nothing)
        replayable === false && throw(
            ArgumentError(
                "Manifest is descriptive only and is not replayable; use a model-level construct_with_manifest method to export a replay recipe."
            ),
        )
        throw(
            ArgumentError(
                "Manifest does not contain a top-level \"recipe\" section required for reconstruction."
            ),
        )
    end

    recipe = manifest["recipe"]
    recipe isa AbstractDict || throw(
        ArgumentError("Manifest recipe must be a dictionary, got $(typeof(recipe)).")
    )

    recipe_type = get(recipe, "type", nothing)
    recipe_type isa AbstractString || throw(
        ArgumentError("Manifest recipe is missing a string \"type\" entry.")
    )
    recipe_type == MODEL_CONSTRUCTOR_RECIPE_TYPE || throw(
        ArgumentError(
            "Unsupported manifest recipe type $(repr(recipe_type)); expected $(repr(MODEL_CONSTRUCTOR_RECIPE_TYPE))."
        ),
    )

    family = get(recipe, "family", nothing)
    family isa AbstractString || throw(
        ArgumentError("Manifest recipe is missing a string \"family\" entry.")
    )
    _model_constructor_family(recipe)

    kwargs_dict = get(recipe, "kwargs", nothing)
    kwargs_dict isa AbstractDict || throw(
        ArgumentError("Manifest recipe is missing a dictionary \"kwargs\" entry.")
    )

    return recipe
end


function _model_constructor_family(recipe::AbstractDict)
    family = recipe["family"]
    family in MODEL_CONSTRUCTOR_FAMILIES && return family
    throw(ArgumentError("Unsupported model constructor family $(repr(family))."))
end

function construct_from_manifest(path::AbstractString; grid=nothing, arch=nothing)
    return construct_from_manifest(JSON.parsefile(path); grid=grid, arch=arch)
end

function _recipe_constructor_kwargs(kwargs_dict::AbstractDict, manifest::AbstractDict; grid=nothing, arch=nothing)
    pairs = Pair{Symbol,Any}[]

    for key in ("phyto_size_structure", "zoo_size_structure", "open_bottom")
        haskey(kwargs_dict, key) && push!(pairs, Symbol(key) => kwargs_dict[key])
    end

    if haskey(kwargs_dict, "parameters")
        parameters = kwargs_dict["parameters"]
        parameters isa AbstractDict || throw(
            ArgumentError(
                "Manifest recipe \"parameters\" entry must be a dictionary, got $(typeof(parameters))."
            ),
        )
        push!(pairs, :parameters => _parameter_kwargs(parameters, manifest))
    end

    if haskey(kwargs_dict, "sinking_tracers")
        push!(pairs, :sinking_tracers => _sinking_tracers_kwargs(kwargs_dict["sinking_tracers"]))
    end

    if haskey(kwargs_dict, "scalar_type")
        push!(pairs, :scalar_type => _decode_manifest_scalar_type(kwargs_dict["scalar_type"]))
    end

    !isnothing(grid) && push!(pairs, :grid => grid)
    !isnothing(arch) && push!(pairs, :arch => arch)

    return (; pairs...)
end


function _sinking_tracers_kwargs(sinking)
    isnothing(sinking) && return nothing

    if sinking isa AbstractDict
        return _dict_to_namedtuple(sinking)
    elseif sinking isa AbstractVector
        pairs = Pair{Symbol,Any}[]
        for item in sinking
            item isa AbstractDict || throw(
                ArgumentError(
                    "Manifest recipe ordered sinking tracer entries must be dictionaries, got $(typeof(item))."
                ),
            )
            haskey(item, "name") && haskey(item, "value") || throw(
                ArgumentError(
                    "Manifest recipe ordered sinking tracer entries must contain \"name\" and \"value\"."
                ),
            )
            push!(pairs, Symbol(item["name"]) => _manifest_value_to_julia(item["value"]))
        end
        return (; pairs...)
    else
        throw(
            ArgumentError(
                "Manifest recipe \"sinking_tracers\" entry must be an ordered vector, dictionary, or nothing, got $(typeof(sinking))."
            ),
        )
    end
end

function _parameter_kwargs(parameters::AbstractDict, manifest::AbstractDict)
    records = get(get(manifest, "resolved", Dict{String,Any}()), "parameters", Dict{String,Any}())
    return (; (
        Symbol(k) => _manifest_parameter_value_to_julia(k, v, records) for
        (k, v) in pairs(parameters)
    )...)
end

function _manifest_parameter_value_to_julia(key, value, records)
    record = records isa AbstractDict ? get(records, string(key), nothing) : nothing
    if record isa AbstractDict && get(record, "shape", nothing) == "matrix"
        return _manifest_matrix(value)
    end
    return _manifest_value_to_julia(value)
end

function _manifest_matrix(value)
    value isa AbstractVector || throw(
        ArgumentError(
            "Manifest matrix parameter must be encoded as a vector of rows, got $(typeof(value))."
        ),
    )
    rows = Any[value[i] for i in eachindex(value)]
    isempty(rows) && return Matrix{Any}(undef, 0, 0)
    all(row -> row isa AbstractVector, rows) || throw(
        ArgumentError("Manifest matrix parameter rows must be vectors."),
    )
    ncols = length(rows[1])
    all(row -> length(row) == ncols, rows) || throw(
        ArgumentError("Manifest matrix parameter rows must all have the same length."),
    )
    return [_manifest_value_to_julia(rows[i][j]) for i in eachindex(rows), j in 1:ncols]
end

function _dict_to_namedtuple(dict::AbstractDict)
    return (; (Symbol(k) => _manifest_value_to_julia(v) for (k, v) in pairs(dict))...)
end

_manifest_value_to_julia(x) = x

function _decode_manifest_scalar_type(x)
    isnothing(x) && return nothing
    x isa Type && x <: Real && return x
    x isa AbstractString || throw(
        ArgumentError(
            "Manifest recipe scalar_type must be a string or concrete Real type, got $(typeof(x))."
        ),
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
function default_model_manifest_context(family::Symbol, resolved; recipe=nothing)
    family_name = string(family)
    return (
        model=(
            id=string(family_name, "/default"),
            family=family_name,
            citation="default",
            tag="default",
        ),
        recipe=recipe,
        resolved=resolved,
    )
end


function default_model_recipe(family::Symbol, resolved)
    family_name = string(family)
    diameters = resolved.plankton_diameters_by_group

    return Dict{String,Any}(
        "type" => MODEL_CONSTRUCTOR_RECIPE_TYPE,
        "family" => family_name,
        "kwargs" => Dict{String,Any}(
            "phyto_size_structure" => diameters["P"],
            "zoo_size_structure" => diameters["Z"],
            "parameters" => resolved.parameter_values,
            "sinking_tracers" => resolved.sinking_tracers_recipe,
            "open_bottom" => resolved.sinking["open_bottom"],
            "scalar_type" => resolved.scalar_type,
        ),
    )
end

function finalize_construction_manifest(context)
    model = context.model
    resolved = context.resolved

    manifest = Dict{String,Any}(
        "schema" => CONSTRUCTION_MANIFEST_SCHEMA,
        "created_at" => string(now(UTC)),
        "agate" => Dict{String,Any}(
            "version" => agate_version(), "julia_version" => string(VERSION)
        ),
        "model" => Dict{String,Any}(
            "id" => model.id,
            "family" => model.family,
            "citation" => model.citation,
            "tag" => model.tag,
        ),
        "resolved" => Dict{String,Any}(
            "tracers" => resolved.tracers,
            "auxiliary_fields" => resolved.auxiliary_fields,
            "parameters" => resolved.parameters,
            "plankton_diameters" => resolved.plankton_diameters,
            "plankton_diameters_by_group" => resolved.plankton_diameters_by_group,
            "scalar_type" => resolved.scalar_type,
            "architecture" => resolved.architecture,
            "has_sinking_velocities" => resolved.has_sinking_velocities,
            "sinking" => resolved.sinking,
        ),
    )

    replayable = hasproperty(context, :recipe) && !isnothing(context.recipe)
    manifest["replayable"] = replayable

    if replayable
        manifest["recipe"] = context.recipe
    end

    return manifest
end
