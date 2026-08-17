using JSON

using ..Factories:
    ParameterDefinition,
    ParameterSpec,
    ConstDefault,
    NoDefault,
    FillDefault,
    DiameterIndexedVectorDefault,
    DiameterIndexedMaterialization

using ..Configuration:
    PFTSpecification,
    DiameterListSpecification,
    DiameterRangeSpecification,
    MatrixDefinition,
    matrix_deriver_identifier,
    matrix_deriver_from_identifier

using ..Library.Allometry:
    ConstantParam,
    AllometricParam,
    allometric_relationship_identifier,
    allometric_relationship_from_identifier

const MODEL_RECIPE_SCHEMA = "agate.model_recipe.v1"
const _RECIPE_DOCUMENT_KEYS = ("schema", "model", "recipe")
const _RECIPE_MODEL_KEYS = ("family",)
const _RECIPE_KEYS = (
    "community",
    "parameter_definitions",
    "parameter_overrides",
    "interaction_definitions",
    "interaction_overrides",
    "ecological_roles",
    "interaction_roles",
    "parameter_roles",
    "auxiliary_fields",
    "sinking_tracers",
    "open_bottom",
    "scalar_type",
)
const _SUPPORTED_SPLITTING = (:linear_splitting, :log_splitting)

function _check_keys(x, allowed, path)
    x isa AbstractDict || throw(ArgumentError("$path must be an object."))
    for key in keys(x)
        key in allowed || throw(ArgumentError("$path has unsupported field $(repr(key))."))
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

function _symbol(x, path)
    return Symbol(_string(x, path))
end

_float_type_id(::Type{Float32}) = "Float32"
_float_type_id(::Type{Float64}) = "Float64"
function _float_type_id(T::Type{<:Real})
    throw(ArgumentError("Recipe serialization supports Float32 and Float64; got $T."))
end

function _decode_float_type(x, path)
    name = _string(x, path)
    name == "Float32" && return Float32
    name == "Float64" && return Float64
    throw(ArgumentError("$path has unsupported floating-point type $(repr(name))."))
end

function _finite_float(x::AbstractFloat)
    isfinite(x) || throw(
        ArgumentError("Recipe serialization supports finite floating-point values; got $x.")
    )
    return Float64(x)
end


function _decoded_array(values, path)
    decoded = Any[_decode_value(v, "$path[$i]") for (i, v) in pairs(values)]
    isempty(decoded) && return Float64[]

    T = promote_type(map(typeof, decoded)...)
    T in (Float64, Int, Bool, Symbol) || throw(
        ArgumentError("$path has unsupported array element types.")
    )
    return T[decoded...]
end

function _decode_matrix(rows, path)
    rows isa AbstractVector || throw(ArgumentError("$path must be an array."))
    isempty(rows) && return Matrix{Float64}(undef, 0, 0)
    all(row -> row isa AbstractVector, rows) || throw(
        ArgumentError("$path must contain arrays.")
    )

    ncols = length(first(rows))
    all(row -> length(row) == ncols, rows) || throw(
        ArgumentError("$path must be rectangular.")
    )
    decoded_rows = [_decoded_array(row, "$path[$i]") for (i, row) in pairs(rows)]
    T = promote_type(map(eltype, decoded_rows)...)
    T in (Float64, Int, Bool, Symbol) || throw(
        ArgumentError("$path has unsupported matrix element types.")
    )
    out = Matrix{T}(undef, length(rows), ncols)
    for i in eachindex(decoded_rows), j in 1:ncols
        out[i, j] = decoded_rows[i][j]
    end
    return out
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

_encode_value(x::Nothing) = nothing
_encode_value(x::Bool) = x
function _encode_value(x::Integer)
    x isa Int || throw(ArgumentError("Recipe serialization supports Int integer values; got $(typeof(x))."))
    return x
end
_encode_value(x::AbstractString) = String(x)
_encode_value(x::Symbol) = Dict{String,Any}("type" => "symbol", "value" => String(x))
_encode_value(x::AbstractFloat) = _finite_float(x)

function _encode_value(x::NamedTuple)
    entries = Any[
        Dict{String,Any}("name" => String(k), "value" => _encode_value(v))
        for (k, v) in pairs(x)
    ]
    return Dict{String,Any}("type" => "named_tuple", "entries" => entries)
end

function _encode_value(x::Tuple)
    return Dict{String,Any}(
        "type" => "tuple", "items" => Any[_encode_value(v) for v in x]
    )
end

function _encode_value(x::AbstractVector)
    x isa Vector || throw(ArgumentError("Recipe serialization supports Vector inputs; got $(typeof(x))."))
    return Any[_encode_value(v) for v in x]
end

function _encode_value(x::AbstractMatrix)
    x isa Matrix || throw(ArgumentError("Recipe serialization supports Matrix inputs; got $(typeof(x))."))
    return Any[
        Any[_encode_value(x[i, j]) for j in axes(x, 2)] for i in axes(x, 1)
    ]
end

function _encode_value(x::PFTSpecification)
    return Dict{String,Any}("type" => "pft_specification", "data" => _encode_value(x.data))
end

function _encode_value(x::DiameterListSpecification)
    return Dict{String,Any}(
        "type" => "diameter_list", "diameters" => _encode_value(x.diameters)
    )
end

function _encode_value(x::DiameterRangeSpecification)
    return Dict{String,Any}(
        "type" => "diameter_range",
        "n" => x.n,
        "min_diameter" => _encode_value(x.min_diameter),
        "max_diameter" => _encode_value(x.max_diameter),
        "splitting" => String(x.splitting),
    )
end

function _encode_value(x::ConstantParam)
    return Dict{String,Any}("law" => "constant", "value" => _encode_value(x.value))
end

function _encode_value(x::AllometricParam)
    return Dict{String,Any}(
        "law" => _identifier_string(
            x.model, allometric_relationship_identifier, "allometric_relationship_identifier"
        ),
        "coefficients" => _encode_value(x.coeffs),
    )
end

function _encode_value(x)
    throw(ArgumentError("Cannot serialize recipe value of type $(typeof(x))."))
end

function _decode_value(x, path)
    x === nothing && return nothing
    x isa Bool && return x
    x isa Int && return x
    x isa AbstractFloat && return _finite_float(x)
    x isa AbstractString && return String(x)
    if x isa AbstractVector
        isempty(x) && return Float64[]
        rows = map(v -> v isa AbstractVector, x)
        all(rows) && return _decode_matrix(x, path)
        any(rows) && throw(ArgumentError("$path cannot mix scalar values and array rows."))
        return _decoded_array(x, path)
    end
    x isa AbstractDict || throw(ArgumentError("$path has unsupported JSON value of type $(typeof(x))."))

    if haskey(x, "law")
        law = _symbol(_required(x, "law", path), "$path.law")
        if law === :constant
            _check_keys(x, ("law", "value"), path)
            return ConstantParam(_decode_value(_required(x, "value", path), "$path.value"))
        end
        _check_keys(x, ("law", "coefficients"), path)
        model = _decode_identifier(
            law, allometric_relationship_from_identifier, "$path.law", "allometric relationship"
        )
        coeffs = _decode_value(_required(x, "coefficients", path), "$path.coefficients")
        coeffs isa NamedTuple || throw(ArgumentError("$path.coefficients must decode to a NamedTuple."))
        return AllometricParam(model, coeffs)
    end

    type_name = _string(_required(x, "type", path), "$path.type")
    if type_name == "symbol"
        _check_keys(x, ("type", "value"), path)
        return _symbol(_required(x, "value", path), "$path.value")
    elseif type_name == "named_tuple"
        _check_keys(x, ("type", "entries"), path)
        return _decode_named_tuple(_required(x, "entries", path), "$path.entries")
    elseif type_name == "tuple"
        _check_keys(x, ("type", "items"), path)
        items = _required(x, "items", path)
        items isa AbstractVector || throw(ArgumentError("$path.items must be an array."))
        return Tuple(_decode_value(v, "$path.items[$i]") for (i, v) in pairs(items))
    elseif type_name == "pft_specification"
        _check_keys(x, ("type", "data"), path)
        return PFTSpecification(_decode_value(_required(x, "data", path), "$path.data"))
    elseif type_name == "diameter_list"
        _check_keys(x, ("type", "diameters"), path)
        diameters = _decode_value(_required(x, "diameters", path), "$path.diameters")
        diameters isa AbstractVector || throw(ArgumentError("$path.diameters must decode to a vector."))
        return DiameterListSpecification(diameters)
    elseif type_name == "diameter_range"
        _check_keys(x, ("type", "n", "min_diameter", "max_diameter", "splitting"), path)
        splitting = _symbol(_required(x, "splitting", path), "$path.splitting")
        splitting in _SUPPORTED_SPLITTING || throw(
            ArgumentError("$path.splitting has unsupported splitting method $(repr(splitting)).")
        )
        return DiameterRangeSpecification(
            _count(_required(x, "n", path), "$path.n"),
            _decode_value(_required(x, "min_diameter", path), "$path.min_diameter"),
            _decode_value(_required(x, "max_diameter", path), "$path.max_diameter"),
            splitting,
        )
    end

    throw(ArgumentError("$path has unsupported semantic type $(repr(type_name))."))
end

function _decode_named_tuple(entries, path)
    entries isa AbstractVector || throw(ArgumentError("$path must be an array."))
    names = Symbol[]
    values = Any[]
    for (i, entry) in pairs(entries)
        entry_path = "$path[$i]"
        entry = _check_keys(entry, ("name", "value"), entry_path)
        name = _symbol(_required(entry, "name", entry_path), "$entry_path.name")
        name in names && throw(ArgumentError("$path contains duplicate name $(repr(name))."))
        push!(names, name)
        push!(values, _decode_value(_required(entry, "value", entry_path), "$entry_path.value"))
    end
    return NamedTuple{Tuple(names)}(Tuple(values))
end

_encode_axes(::Nothing) = nothing
_encode_axes(axis::Symbol) = String(axis)
_encode_axes(axes::Tuple) = String[String(axis) for axis in axes]

function _decode_axes(x, path)
    x === nothing && return nothing
    x isa AbstractString && return Symbol(x)
    x isa AbstractVector && length(x) == 2 || throw(
        ArgumentError("$path must be null, a string, or a two-element array.")
    )
    return (_symbol(x[1], "$path[1]"), _symbol(x[2], "$path[2]"))
end

_encode_materialization(::Nothing) = nothing
function _encode_materialization(materialization::DiameterIndexedMaterialization)
    return Dict{String,Any}(
        "role" => isnothing(materialization.role) ? nothing : String(materialization.role),
        "fill_value" => _encode_value(materialization.fill_value),
    )
end

function _decode_materialization(x, path)
    x === nothing && return nothing
    x = _check_keys(x, ("role", "fill_value"), path)
    role = _required(x, "role", path)
    return DiameterIndexedMaterialization(
        isnothing(role) ? nothing : _symbol(role, "$path.role");
        fill_value=_decode_value(_required(x, "fill_value", path), "$path.fill_value"),
    )
end

_encode_default(::NoDefault) = nothing
function _encode_default(provider::Union{ConstDefault,FillDefault})
    return Dict{String,Any}("value" => _encode_value(provider.value))
end
function _encode_default(provider::DiameterIndexedVectorDefault)
    return Dict{String,Any}(
        "value" => _encode_value(provider.value),
        "role" => String(provider.role),
        "fill_value" => _encode_value(provider.default),
    )
end

function _decode_default(x, shape::Symbol, path)
    x === nothing && return NoDefault()
    x isa AbstractDict || throw(ArgumentError("$path must be an object or null."))
    if haskey(x, "role")
        _check_keys(x, ("value", "role", "fill_value"), path)
        shape === :vector || throw(ArgumentError("$path diameter-indexed defaults require vector shape."))
        return DiameterIndexedVectorDefault(
            _decode_value(_required(x, "value", path), "$path.value"),
            _symbol(_required(x, "role", path), "$path.role");
            default=_decode_value(_required(x, "fill_value", path), "$path.fill_value"),
        )
    end
    _check_keys(x, ("value",), path)
    value = _decode_value(_required(x, "value", path), "$path.value")
    return shape === :scalar ? ConstDefault(value) : FillDefault(value)
end

function _encode_parameter_definitions(definitions)
    return Any[
        Dict{String,Any}(
            "name" => String(definition.spec.name),
            "shape" => String(definition.spec.shape),
            "axes" => _encode_axes(definition.spec.axes),
            "materialization" => _encode_materialization(definition.spec.materialization),
            "default" => _encode_default(definition.default),
        ) for definition in definitions
    ]
end

function _decode_parameter_definitions(x, path)
    x isa AbstractVector || throw(ArgumentError("$path must be an array."))
    return Tuple(begin
        item_path = "$path[$i]"
        item = _check_keys(
            value, ("name", "shape", "axes", "materialization", "default"), item_path
        )
        shape = _symbol(_required(item, "shape", item_path), "$item_path.shape")
        shape in (:scalar, :vector, :matrix) || throw(
            ArgumentError("$item_path.shape has unsupported shape $(repr(shape)).")
        )
        spec = ParameterSpec(
            _symbol(_required(item, "name", item_path), "$item_path.name"),
            shape;
            axes=_decode_axes(_required(item, "axes", item_path), "$item_path.axes"),
            materialization=_decode_materialization(
                _required(item, "materialization", item_path), "$item_path.materialization"
            ),
        )
        ParameterDefinition(
            spec, _decode_default(_required(item, "default", item_path), shape, "$item_path.default")
        )
    end for (i, value) in pairs(x))
end

function _encode_interaction_definitions(definitions::NamedTuple)
    return Any[
        Dict{String,Any}(
            "name" => String(name),
            "deriver" => _identifier_string(
                definition.deriver, matrix_deriver_identifier, "matrix_deriver_identifier"
            ),
            "deps" => String[String(dep) for dep in definition.deps],
        ) for (name, definition) in pairs(definitions)
    ]
end

function _decode_interaction_definitions(x, path)
    x isa AbstractVector || throw(ArgumentError("$path must be an array."))
    names = Symbol[]
    values = Any[]
    for (i, item) in pairs(x)
        item_path = "$path[$i]"
        item = _check_keys(item, ("name", "deriver", "deps"), item_path)
        name = _symbol(_required(item, "name", item_path), "$item_path.name")
        name in names && throw(ArgumentError("$path contains duplicate name $(repr(name))."))
        deriver_name = _symbol(_required(item, "deriver", item_path), "$item_path.deriver")
        deriver = _decode_identifier(
            deriver_name,
            matrix_deriver_from_identifier,
            "$item_path.deriver",
            "interaction-matrix deriver",
        )
        deps_value = _required(item, "deps", item_path)
        deps_value isa AbstractVector || throw(ArgumentError("$item_path.deps must be an array."))
        deps = Tuple(_symbol(dep, "$item_path.deps[$j]") for (j, dep) in pairs(deps_value))
        push!(names, name)
        push!(values, MatrixDefinition(deriver; deps))
    end
    return NamedTuple{Tuple(names)}(Tuple(values))
end

function _encode_community(community::NamedTuple)
    return Any[
        Dict{String,Any}("name" => String(name), "spec" => _encode_value(spec))
        for (name, spec) in pairs(community)
    ]
end

function _decode_community(x, path)
    x isa AbstractVector || throw(ArgumentError("$path must be an array."))
    names = Symbol[]
    specs = Any[]
    for (i, item) in pairs(x)
        item_path = "$path[$i]"
        item = _check_keys(item, ("name", "spec"), item_path)
        name = _symbol(_required(item, "name", item_path), "$item_path.name")
        name in names && throw(ArgumentError("$path contains duplicate group $(repr(name))."))
        spec = _decode_value(_required(item, "spec", item_path), "$item_path.spec")
        spec isa NamedTuple || throw(ArgumentError("$item_path.spec must decode to a NamedTuple."))
        push!(names, name)
        push!(specs, spec)
    end
    return NamedTuple{Tuple(names)}(Tuple(specs))
end

"""Encode a `ModelRecipe` as a JSON-compatible recipe document."""
function encode_recipe(recipe::ModelRecipe)
    data = Dict{String,Any}(
        "community" => _encode_community(recipe.community),
        "parameter_definitions" => _encode_parameter_definitions(recipe.parameter_definitions),
        "parameter_overrides" => _encode_value(recipe.parameter_overrides),
        "interaction_definitions" => _encode_interaction_definitions(recipe.interaction_definitions),
        "interaction_overrides" => _encode_value(recipe.interaction_overrides),
        "ecological_roles" => _encode_value(recipe.ecological_roles),
        "interaction_roles" => _encode_value(recipe.interaction_roles),
        "parameter_roles" => _encode_value(recipe.parameter_roles),
        "auxiliary_fields" => _encode_value(recipe.auxiliary_fields),
        "sinking_tracers" => isnothing(recipe.sinking_tracers) ? nothing : _encode_value(recipe.sinking_tracers),
        "open_bottom" => recipe.open_bottom,
        "scalar_type" => _float_type_id(recipe.scalar_type),
    )

    return Dict{String,Any}(
        "schema" => MODEL_RECIPE_SCHEMA,
        "model" => Dict{String,Any}("family" => String(recipe.family)),
        "recipe" => data,
    )
end

"""Decode a recipe document into a `ModelRecipe`."""
function decode_recipe(document::AbstractDict)
    document = _check_keys(document, _RECIPE_DOCUMENT_KEYS, "Recipe document")
    schema = _string(_required(document, "schema", "Recipe document"), "Recipe document.schema")
    schema == MODEL_RECIPE_SCHEMA || throw(
        ArgumentError("Unsupported Agate recipe schema $(repr(schema)); expected $(repr(MODEL_RECIPE_SCHEMA)).")
    )

    model = _check_keys(_required(document, "model", "Recipe document"), _RECIPE_MODEL_KEYS, "Recipe document.model")
    family = _symbol(_required(model, "family", "Recipe document.model"), "Recipe document.model.family")
    recipe_factory(Val(family))

    recipe = _check_keys(_required(document, "recipe", "Recipe document"), _RECIPE_KEYS, "Recipe document.recipe")
    for key in _RECIPE_KEYS
        _required(recipe, key, "Recipe document.recipe")
    end

    community = _decode_community(recipe["community"], "Recipe document.recipe.community")
    parameter_definitions = _decode_parameter_definitions(recipe["parameter_definitions"], "Recipe document.recipe.parameter_definitions")
    parameter_overrides = _decode_value(recipe["parameter_overrides"], "Recipe document.recipe.parameter_overrides")
    parameter_overrides isa NamedTuple || throw(ArgumentError("Recipe document.recipe.parameter_overrides must decode to a NamedTuple."))
    interaction_definitions = _decode_interaction_definitions(recipe["interaction_definitions"], "Recipe document.recipe.interaction_definitions")
    interaction_overrides = _decode_value(recipe["interaction_overrides"], "Recipe document.recipe.interaction_overrides")
    interaction_overrides isa NamedTuple || throw(ArgumentError("Recipe document.recipe.interaction_overrides must decode to a NamedTuple."))
    all(value -> value isa AbstractMatrix, values(interaction_overrides)) || throw(
        ArgumentError("Recipe document.recipe.interaction_overrides values must be matrices.")
    )
    ecological_roles = _decode_value(recipe["ecological_roles"], "Recipe document.recipe.ecological_roles")
    interaction_roles = _decode_value(recipe["interaction_roles"], "Recipe document.recipe.interaction_roles")
    parameter_roles = _decode_value(recipe["parameter_roles"], "Recipe document.recipe.parameter_roles")
    auxiliary_fields = _decode_value(recipe["auxiliary_fields"], "Recipe document.recipe.auxiliary_fields")
    sinking_tracers = isnothing(recipe["sinking_tracers"]) ? nothing : _decode_value(recipe["sinking_tracers"], "Recipe document.recipe.sinking_tracers")

    ecological_roles isa NamedTuple || throw(ArgumentError("Recipe document.recipe.ecological_roles must decode to a NamedTuple."))
    interaction_roles isa NamedTuple || throw(ArgumentError("Recipe document.recipe.interaction_roles must decode to a NamedTuple."))
    parameter_roles isa NamedTuple || throw(ArgumentError("Recipe document.recipe.parameter_roles must decode to a NamedTuple."))
    auxiliary_fields isa Tuple || throw(ArgumentError("Recipe document.recipe.auxiliary_fields must decode to a Tuple."))
    !isnothing(sinking_tracers) && !(sinking_tracers isa NamedTuple) && throw(ArgumentError("Recipe document.recipe.sinking_tracers must decode to a NamedTuple or null."))

    return ModelRecipe(
        family,
        community,
        parameter_definitions,
        parameter_overrides,
        interaction_definitions,
        interaction_overrides,
        ecological_roles,
        interaction_roles,
        parameter_roles,
        auxiliary_fields,
        sinking_tracers,
        _boolean(recipe["open_bottom"], "Recipe document.recipe.open_bottom"),
        _decode_float_type(recipe["scalar_type"], "Recipe document.recipe.scalar_type"),
    )
end

function decode_recipe(document)
    throw(ArgumentError("Expected an Agate recipe dictionary, got $(typeof(document))."))
end

"""Write a recipe document to `path` as pretty-printed JSON."""
function export_recipe(path::AbstractString, recipe::ModelRecipe)
    document = encode_recipe(recipe)
    open(path, "w") do io
        JSON.print(io, document, 4)
        println(io)
    end
    return path
end

"""Read and decode a recipe document from `path`."""
import_recipe(path::AbstractString) = decode_recipe(JSON.parsefile(path))
