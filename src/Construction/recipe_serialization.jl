using JSON

using ..Configuration: DiameterListSpecification, DiameterRangeSpecification
using ..Library.Allometry:
    ConstantParam,
    AllometricParam,
    allometric_relationship_identifier,
    allometric_relationship_from_identifier

const MODEL_RECIPE_SCHEMA = "agate.model_recipe.v1"
const _RECIPE_DOCUMENT_KEYS = (
    "schema", "family", "definition_version", "realization", "provenance", "content_hash"
)
const _REALIZATION_KEYS = (
    "plankton_pfts", "pft_size_structures", "parameter_overrides", "sinking_tracers", "open_bottom"
)
const _SUPPORTED_SPACING = (:linear, :log)

function _check_keys(x, allowed, path)
    x isa AbstractDict || throw(ArgumentError("$path must be an object."))
    for key in keys(x)
        key in allowed || throw(ArgumentError("$path has unsupported field $(repr(key))."))
    end
    return x
end

function _complete_object(x, keys, path)
    x = _check_keys(x, keys, path)
    for key in keys
        _required(x, key, path)
    end
    return x
end

function _required(x::AbstractDict, key, path)
    haskey(x, key) || throw(ArgumentError("$path is missing required field $(repr(key))."))
    return x[key]
end

function _string(x, path)
    x isa AbstractString || throw(ArgumentError("$path must be a string."))
    return String(x)
end

function _boolean(x, path)
    x isa Bool || throw(ArgumentError("$path must be a boolean."))
    return x
end

function _count(x, path)
    x isa Integer && !(x isa Bool) || throw(ArgumentError("$path must be an integer."))
    x > 0 || throw(ArgumentError("$path must be a positive integer."))
    return Int(x)
end

_symbol(x, path) = Symbol(_string(x, path))

function _version(x, path)
    text = _string(x, path)
    try
        return VersionNumber(text)
    catch
        throw(ArgumentError("$path must be a valid version number."))
    end
end

function _finite_float(x::AbstractFloat)
    isfinite(x) || throw(
        ArgumentError("Recipe serialization supports finite floating-point values; got $x.")
    )
    return Float64(x)
end

function _identifier_string(value, identifier_function, function_name)
    identifier = identifier_function(value)
    identifier isa Symbol || throw(
        ArgumentError("$function_name must return a Symbol; got $(typeof(identifier)).")
    )
    return String(identifier)
end

function _decode_identifier(identifier::Symbol, decoder, path, kind)
    try
        return decoder(Val(identifier))
    catch err
        err isa ArgumentError || rethrow()
        throw(ArgumentError("$path has unsupported $kind $(repr(identifier))."))
    end
end

# Narrow codec for authored parameter values. The durable recipe intentionally does not
# serialize arbitrary Julia values or the component/process definition tree.
_encode_parameter_value(x::Nothing) = nothing
_encode_parameter_value(x::Bool) = x
function _encode_parameter_value(x::Integer)
    return Int(x)
end
_encode_parameter_value(x::AbstractFloat) = _finite_float(x)
_encode_parameter_value(x::Symbol) = Dict{String,Any}("kind" => "symbol", "value" => String(x))

function _encode_parameter_value(x::NamedTuple)
    return Dict{String,Any}(
        "kind" => "named",
        "entries" => Any[
            Dict{String,Any}("name" => String(name), "value" => _encode_parameter_value(value))
            for (name, value) in pairs(x)
        ],
    )
end

function _encode_parameter_value(x::AbstractVector)
    all(value -> value isa Real && isfinite(value), x) || throw(
        ArgumentError("Recipe parameter vectors must contain finite numeric values.")
    )
    return Any[_encode_parameter_value(value) for value in x]
end

function _encode_parameter_value(x::AbstractMatrix)
    all(value -> value isa Real && isfinite(value), x) || throw(
        ArgumentError("Recipe parameter matrices must contain finite numeric values.")
    )
    return Any[
        Any[_encode_parameter_value(x[i, j]) for j in axes(x, 2)] for i in axes(x, 1)
    ]
end

function _encode_parameter_value(x::ConstantParam)
    return Dict{String,Any}("law" => "constant", "value" => _encode_parameter_value(x.value))
end

function _encode_parameter_value(x::AllometricParam)
    return Dict{String,Any}(
        "law" => _identifier_string(
            x.model, allometric_relationship_identifier, "allometric_relationship_identifier"
        ),
        "coefficients" => _encode_parameter_value(x.coeffs),
    )
end

function _encode_parameter_value(x)
    throw(ArgumentError("Cannot serialize recipe parameter value of type $(typeof(x))."))
end

function _decoded_parameter_array(values, path)
    decoded = Any[_decode_parameter_value(value, "$path[$i]") for (i, value) in pairs(values)]
    isempty(decoded) && return Float64[]
    rows = map(value -> value isa AbstractVector, decoded)
    any(rows) && !all(rows) && throw(
        ArgumentError("$path cannot mix scalar values and array rows.")
    )
    if all(rows)
        ncols = length(first(decoded))
        all(row -> length(row) == ncols, decoded) ||
            throw(ArgumentError("$path must be rectangular."))
        T = promote_type(map(eltype, decoded)...)
        T in (Float64, Int, Bool, Symbol) || throw(
            ArgumentError("$path has unsupported matrix element types.")
        )
        out = Matrix{T}(undef, length(decoded), ncols)
        for i in eachindex(decoded), j in 1:ncols
            out[i, j] = decoded[i][j]
        end
        return out
    end

    T = promote_type(map(typeof, decoded)...)
    T in (Float64, Int, Bool, Symbol) || throw(
        ArgumentError("$path has unsupported vector element types.")
    )
    return T[decoded...]
end

function _decode_named_parameter_values(entries, path)
    entries isa AbstractVector || throw(ArgumentError("$path must be an array."))
    names = Symbol[]
    values = Any[]
    for (i, entry) in pairs(entries)
        entry_path = "$path[$i]"
        entry = _complete_object(entry, ("name", "value"), entry_path)
        name = _symbol(entry["name"], "$entry_path.name")
        name in names && throw(ArgumentError("$path contains duplicate name $(repr(name))."))
        push!(names, name)
        push!(values, _decode_parameter_value(entry["value"], "$entry_path.value"))
    end
    return NamedTuple{Tuple(names)}(Tuple(values))
end

function _decode_parameter_value(x, path)
    x === nothing && return nothing
    x isa Bool && return x
    x isa Integer && !(x isa Bool) && return Int(x)
    x isa AbstractFloat && return _finite_float(x)
    x isa AbstractVector && return _decoded_parameter_array(x, path)
    x isa AbstractDict || throw(
        ArgumentError("$path has unsupported JSON value of type $(typeof(x)).")
    )

    if haskey(x, "law")
        law = _symbol(_required(x, "law", path), "$path.law")
        if law === :constant
            _complete_object(x, ("law", "value"), path)
            return ConstantParam(_decode_parameter_value(x["value"], "$path.value"))
        end
        _complete_object(x, ("law", "coefficients"), path)
        model = _decode_identifier(
            law, allometric_relationship_from_identifier, "$path.law", "allometric relationship"
        )
        coefficients = _decode_parameter_value(x["coefficients"], "$path.coefficients")
        coefficients isa NamedTuple || throw(
            ArgumentError("$path.coefficients must decode to a NamedTuple.")
        )
        return AllometricParam(model, coefficients)
    end

    kind = _string(_required(x, "kind", path), "$path.kind")
    if kind == "symbol"
        _complete_object(x, ("kind", "value"), path)
        return _symbol(x["value"], "$path.value")
    elseif kind == "named"
        _complete_object(x, ("kind", "entries"), path)
        return _decode_named_parameter_values(x["entries"], "$path.entries")
    end
    throw(ArgumentError("$path has unsupported parameter value kind $(repr(kind))."))
end

function _encode_parameter_overrides(overrides::NamedTuple)
    return Any[
        Dict{String,Any}("name" => String(name), "value" => _encode_parameter_value(value))
        for (name, value) in pairs(overrides)
    ]
end

function _decode_parameter_overrides(x, path)
    return _decode_named_parameter_values(x, path)
end

_encode_diameter_specification(::Nothing) = nothing

function _encode_diameter_specification(spec::DiameterListSpecification)
    return Dict{String,Any}(
        "kind" => "list",
        "diameters" => Any[_finite_float(float(value)) for value in spec.diameters],
    )
end

function _encode_diameter_specification(spec::DiameterRangeSpecification)
    return Dict{String,Any}(
        "kind" => "range",
        "n" => Int(spec.n),
        "min_esd" => _finite_float(float(spec.min_diameter)),
        "max_esd" => _finite_float(float(spec.max_diameter)),
        "spacing" => String(spec.spacing),
    )
end

function _decode_diameter_specification(x, path)
    x === nothing && return nothing
    x isa AbstractDict || throw(ArgumentError("$path must be an object or null."))
    kind = _string(_required(x, "kind", path), "$path.kind")
    if kind == "list"
        _complete_object(x, ("kind", "diameters"), path)
        values = x["diameters"]
        values isa AbstractVector || throw(ArgumentError("$path.diameters must be an array."))
        diameters = Float64[]
        for (i, value) in pairs(values)
            value isa Real && !(value isa Bool) || throw(
                ArgumentError("$path.diameters[$i] must be numeric.")
            )
            push!(diameters, _finite_float(float(value)))
        end
        return DiameterListSpecification(diameters)
    elseif kind == "range"
        _complete_object(x, ("kind", "n", "min_esd", "max_esd", "spacing"), path)
        spacing = _symbol(x["spacing"], "$path.spacing")
        spacing in _SUPPORTED_SPACING || throw(
            ArgumentError("$path.spacing has unsupported spacing method $(repr(spacing)).")
        )
        min_esd = x["min_esd"]
        max_esd = x["max_esd"]
        min_esd isa Real && !(min_esd isa Bool) ||
            throw(ArgumentError("$path.min_esd must be numeric."))
        max_esd isa Real && !(max_esd isa Bool) ||
            throw(ArgumentError("$path.max_esd must be numeric."))
        return DiameterRangeSpecification(
            _count(x["n"], "$path.n"),
            _finite_float(float(min_esd)),
            _finite_float(float(max_esd)),
            spacing,
        )
    end
    throw(ArgumentError("$path has unsupported diameter kind $(repr(kind))."))
end

function _encode_plankton_pfts(plankton_pfts::NamedTuple)
    return Any[
        Dict{String,Any}(
            "plankton" => String(plankton),
            "pfts" => Any[String(pft) for pft in pfts],
        )
        for (plankton, pfts) in pairs(plankton_pfts)
    ]
end

function _decode_plankton_pfts(x, path)
    x isa AbstractVector || throw(ArgumentError("$path must be an array."))
    plankton_names = Symbol[]
    values = Any[]
    for (i, item) in pairs(x)
        item_path = "$path[$i]"
        item = _complete_object(item, ("plankton", "pfts"), item_path)
        plankton = _symbol(item["plankton"], "$item_path.plankton")
        plankton in plankton_names && throw(
            ArgumentError("$path contains duplicate plankton $(repr(plankton)).")
        )
        pfts_raw = item["pfts"]
        pfts_raw isa AbstractVector || throw(ArgumentError("$item_path.pfts must be an array."))
        isempty(pfts_raw) && throw(ArgumentError("$item_path.pfts cannot be empty."))
        pfts = Tuple(_symbol(pft, "$item_path.pfts[$j]") for (j, pft) in pairs(pfts_raw))
        length(unique(pfts)) == length(pfts) || throw(
            ArgumentError("$item_path.pfts contains duplicate PFT names.")
        )
        push!(plankton_names, plankton)
        push!(values, pfts)
    end
    return NamedTuple{Tuple(plankton_names)}(Tuple(values))
end

function _encode_pft_size_structures(pft_size_structures::NamedTuple)
    return Any[
        Dict{String,Any}(
            "pft" => String(pft),
            "diameters" => _encode_diameter_specification(specification),
        )
        for (pft, specification) in pairs(pft_size_structures)
    ]
end

function _decode_pft_size_structures(x, path)
    x isa AbstractVector || throw(ArgumentError("$path must be an array."))
    pfts = Symbol[]
    specifications = Any[]
    for (i, item) in pairs(x)
        item_path = "$path[$i]"
        item = _complete_object(item, ("pft", "diameters"), item_path)
        pft = _symbol(item["pft"], "$item_path.pft")
        pft in pfts && throw(ArgumentError("$path contains duplicate PFT $(repr(pft))."))
        push!(pfts, pft)
        push!(
            specifications,
            _decode_diameter_specification(item["diameters"], "$item_path.diameters"),
        )
    end
    return NamedTuple{Tuple(pfts)}(Tuple(specifications))
end

function _validate_realization(plankton_pfts::NamedTuple, pft_size_structures::NamedTuple)
    referenced = Symbol[pft for pfts in values(plankton_pfts) for pft in pfts]
    length(unique(referenced)) == length(referenced) || throw(
        ArgumentError("Recipe plankton realization assigns one PFT more than once.")
    )
    Set(referenced) == Set(keys(pft_size_structures)) || throw(
        ArgumentError("Recipe PFT size structures must match the PFTs assigned to plankton exactly.")
    )
    return nothing
end

function _encode_realization(recipe::ModelRecipe)
    _validate_realization(recipe.plankton_pfts, recipe.pft_size_structures)
    return Dict{String,Any}(
        "plankton_pfts" => _encode_plankton_pfts(recipe.plankton_pfts),
        "pft_size_structures" => _encode_pft_size_structures(recipe.pft_size_structures),
        "parameter_overrides" => _encode_parameter_overrides(recipe.parameter_overrides),
        "sinking_tracers" => isnothing(recipe.sinking_tracers) ? nothing :
                             _encode_parameter_overrides(recipe.sinking_tracers),
        "open_bottom" => recipe.open_bottom,
    )
end

function _decode_realization(x, path)
    realization = _complete_object(x, _REALIZATION_KEYS, path)
    plankton_pfts = _decode_plankton_pfts(
        realization["plankton_pfts"], "$path.plankton_pfts"
    )
    pft_size_structures = _decode_pft_size_structures(realization["pft_size_structures"], "$path.pft_size_structures")
    _validate_realization(plankton_pfts, pft_size_structures)
    parameter_overrides = _decode_parameter_overrides(
        realization["parameter_overrides"], "$path.parameter_overrides"
    )
    sinking_tracers = isnothing(realization["sinking_tracers"]) ? nothing :
                      _decode_parameter_overrides(
                          realization["sinking_tracers"], "$path.sinking_tracers"
                      )
    open_bottom = _boolean(realization["open_bottom"], "$path.open_bottom")
    return (; plankton_pfts, pft_size_structures, parameter_overrides, sinking_tracers, open_bottom)
end

"""Encode a versioned family recipe with a scientific content hash and package provenance."""
function encode_recipe(recipe::ModelRecipe)
    realization = _encode_realization(recipe)
    return Dict{String,Any}(
        "schema" => MODEL_RECIPE_SCHEMA,
        "family" => String(recipe.family),
        "definition_version" => string(recipe.definition_version),
        "realization" => realization,
        "provenance" => _recipe_provenance(recipe),
        "content_hash" => _recipe_hash(recipe.family, recipe.definition_version, realization),
    )
end

"""Decode a recipe document, verifying its hash, family registration, and definition version."""
function decode_recipe(document::AbstractDict)
    document = _complete_object(document, _RECIPE_DOCUMENT_KEYS, "Recipe document")
    schema = _string(document["schema"], "Recipe document.schema")
    schema == MODEL_RECIPE_SCHEMA || throw(
        ArgumentError(
            "Unsupported Agate recipe schema $(repr(schema)); supported schema is " *
            "$(repr(MODEL_RECIPE_SCHEMA))."
        )
    )

    family_id_value = _symbol(document["family"], "Recipe document.family")
    version = _version(document["definition_version"], "Recipe document.definition_version")
    realization_data = _complete_object(
        document["realization"], _REALIZATION_KEYS, "Recipe document.realization"
    )
    recorded_hash = _string(document["content_hash"], "Recipe document.content_hash")
    recorded_hash == _recipe_hash(family_id_value, version, realization_data) || throw(
        ArgumentError("Recipe document.content_hash does not match the serialized recipe content.")
    )

    family = _resolve_recipe_family(family_id_value, version)

    realization = _decode_realization(realization_data, "Recipe document.realization")
    plankton_pfts, pft_size_structures = _canonical_recipe_realization(
        family, realization.plankton_pfts, realization.pft_size_structures
    )
    decoded = ModelRecipe(
        family_id_value,
        version,
        plankton_pfts,
        pft_size_structures,
        realization.parameter_overrides,
        realization.sinking_tracers,
        realization.open_bottom,
    )
    provenance = _decode_provenance(document["provenance"], "Recipe document.provenance")
    _check_recipe_provenance(decoded, provenance)
    return decoded
end

function decode_recipe(document)
    throw(ArgumentError("Expected an Agate recipe dictionary, got $(typeof(document))."))
end

"""Write a recipe document to `path` as pretty-printed JSON."""
function export_recipe(path::AbstractString, recipe::ModelRecipe)
    open(path, "w") do io
        JSON.json(io, encode_recipe(recipe); pretty=4, omit_empty=false)
        println(io)
    end
    return path
end

"""Read and decode a recipe document from `path`."""
import_recipe(path::AbstractString) = decode_recipe(JSON.parsefile(path))
