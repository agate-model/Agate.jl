using JSON

using ..Factories:
    AbstractBGCFactory,
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

_scalar_type_id(::Type{Float32}) = "Float32"
_scalar_type_id(::Type{Float64}) = "Float64"
function _scalar_type_id(T::Type{<:Real})
    throw(ArgumentError("Recipe serialization supports scalar types Float32 and Float64; got $T."))
end

function _decode_scalar_type(x, path)
    name = _string(x, path)
    name == "Float32" && return Float32
    name == "Float64" && return Float64
    throw(ArgumentError("$path has unsupported scalar type $(repr(name))."))
end

function _float_format(::Type{Float32})
    return "Float32"
end
function _float_format(::Type{Float64})
    return "Float64"
end
function _float_format(T::Type{<:AbstractFloat})
    throw(ArgumentError("Recipe serialization supports Float32 and Float64 values; got $T."))
end

function _float_payload(x::AbstractFloat)
    isnan(x) && return "NaN"
    x == Inf && return "Inf"
    x == -Inf && return "-Inf"
    return x
end

function _decode_float_payload(x, ::Type{T}, path) where {T<:AbstractFloat}
    if x isa Number && !(x isa Bool)
        return T(x)
    elseif x isa AbstractString
        x == "NaN" && return T(NaN)
        x == "Inf" && return T(Inf)
        x == "-Inf" && return T(-Inf)
    end
    throw(ArgumentError("$path must be a finite number or one of \"NaN\", \"Inf\", \"-Inf\"."))
end

function _element_type_id(T::Type)
    T === Float32 && return "Float32"
    T === Float64 && return "Float64"
    T === Int && return "Int"
    T === Bool && return "Bool"
    T === Symbol && return "Symbol"
    throw(ArgumentError("Recipe serialization does not support array element type $T."))
end

function _decode_element_type(x, path)
    name = _string(x, path)
    name == "Float32" && return Float32
    name == "Float64" && return Float64
    name == "Int" && return Int
    name == "Bool" && return Bool
    name == "Symbol" && return Symbol
    throw(ArgumentError("$path has unsupported array element type $(repr(name))."))
end

_encode_array_item(x::AbstractFloat) = _float_payload(x)
_encode_array_item(x::Symbol) = String(x)
_encode_array_item(x) = _encode_value(x)

function _decode_array_item(x, ::Type{T}, path) where {T<:AbstractFloat}
    return _decode_float_payload(x, T, path)
end
function _decode_array_item(x, ::Type{Int}, path)
    x isa Integer && !(x isa Bool) || throw(ArgumentError("$path must be an integer."))
    return Int(x)
end
function _decode_array_item(x, ::Type{Bool}, path)
    return _boolean(x, path)
end
function _decode_array_item(x, ::Type{Symbol}, path)
    return _symbol(x, path)
end
function _relationship_id(model)
    identifier = allometric_relationship_identifier(model)
    identifier isa Symbol || throw(
        ArgumentError(
            "allometric_relationship_identifier must return a Symbol; got $(typeof(identifier))."
        ),
    )
    return String(identifier)
end

function _decode_relationship(identifier::Symbol, path)
    try
        return allometric_relationship_from_identifier(Val(identifier))
    catch err
        err isa ArgumentError || rethrow()
        throw(ArgumentError("$path has unsupported allometric relationship $(repr(identifier))."))
    end
end

function _deriver_id(deriver)
    identifier = matrix_deriver_identifier(deriver)
    identifier isa Symbol || throw(
        ArgumentError(
            "matrix_deriver_identifier must return a Symbol; got $(typeof(identifier))."
        ),
    )
    return String(identifier)
end

function _decode_deriver(identifier::Symbol, path)
    try
        return matrix_deriver_from_identifier(Val(identifier))
    catch err
        err isa ArgumentError || rethrow()
        throw(ArgumentError("$path has unsupported interaction-matrix deriver $(repr(identifier))."))
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
function _encode_value(x::AbstractFloat)
    return Dict{String,Any}(
        "type" => "float",
        "format" => _float_format(typeof(x)),
        "value" => _float_payload(x),
    )
end

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
    return Dict{String,Any}(
        "type" => "vector",
        "element_type" => _element_type_id(eltype(x)),
        "values" => Any[_encode_array_item(v) for v in x],
    )
end

function _encode_value(x::AbstractMatrix)
    x isa Matrix || throw(ArgumentError("Recipe serialization supports Matrix inputs; got $(typeof(x))."))
    rows = Any[
        Any[_encode_array_item(x[i, j]) for j in axes(x, 2)] for i in axes(x, 1)
    ]
    return Dict{String,Any}(
        "type" => "matrix",
        "element_type" => _element_type_id(eltype(x)),
        "rows" => rows,
    )
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
    return Dict{String,Any}("type" => "constant_param", "value" => _encode_value(x.value))
end

function _encode_value(x::AllometricParam)
    return Dict{String,Any}(
        "type" => "allometric_param",
        "relationship" => _relationship_id(x.model),
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
    x isa AbstractString && return String(x)
    x isa AbstractDict || throw(ArgumentError("$path has unsupported JSON value of type $(typeof(x))."))

    type_name = _string(_required(x, "type", path), "$path.type")
    if type_name == "symbol"
        _check_keys(x, ("type", "value"), path)
        return _symbol(_required(x, "value", path), "$path.value")
    elseif type_name == "float"
        _check_keys(x, ("type", "format", "value"), path)
        T = _decode_scalar_type(_required(x, "format", path), "$path.format")
        return _decode_float_payload(_required(x, "value", path), T, "$path.value")
    elseif type_name == "named_tuple"
        _check_keys(x, ("type", "entries"), path)
        return _decode_named_tuple(_required(x, "entries", path), "$path.entries")
    elseif type_name == "tuple"
        _check_keys(x, ("type", "items"), path)
        items = _required(x, "items", path)
        items isa AbstractVector || throw(ArgumentError("$path.items must be an array."))
        return Tuple(_decode_value(v, "$path.items[$i]") for (i, v) in pairs(items))
    elseif type_name == "vector"
        _check_keys(x, ("type", "element_type", "values"), path)
        T = _decode_element_type(_required(x, "element_type", path), "$path.element_type")
        values = _required(x, "values", path)
        values isa AbstractVector || throw(ArgumentError("$path.values must be an array."))
        return T[_decode_array_item(v, T, "$path.values[$i]") for (i, v) in pairs(values)]
    elseif type_name == "matrix"
        _check_keys(x, ("type", "element_type", "rows"), path)
        T = _decode_element_type(_required(x, "element_type", path), "$path.element_type")
        rows = _required(x, "rows", path)
        rows isa AbstractVector || throw(ArgumentError("$path.rows must be an array."))
        isempty(rows) && return Matrix{T}(undef, 0, 0)
        all(row -> row isa AbstractVector, rows) || throw(ArgumentError("$path.rows must contain arrays."))
        ncols = length(first(rows))
        all(row -> length(row) == ncols, rows) || throw(ArgumentError("$path.rows must be rectangular."))
        out = Matrix{T}(undef, length(rows), ncols)
        for i in eachindex(rows), j in 1:ncols
            out[i, j] = _decode_array_item(rows[i][j], T, "$path.rows[$i][$j]")
        end
        return out
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
    elseif type_name == "constant_param"
        _check_keys(x, ("type", "value"), path)
        return ConstantParam(_decode_value(_required(x, "value", path), "$path.value"))
    elseif type_name == "allometric_param"
        _check_keys(x, ("type", "relationship", "coefficients"), path)
        relationship = _symbol(_required(x, "relationship", path), "$path.relationship")
        model = _decode_relationship(relationship, "$path.relationship")
        coeffs = _decode_value(_required(x, "coefficients", path), "$path.coefficients")
        coeffs isa NamedTuple || throw(ArgumentError("$path.coefficients must decode to a NamedTuple."))
        return AllometricParam(model, coeffs)
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

function _encode_parameter_spec(spec::ParameterSpec)
    axes = if spec.axes === nothing
        nothing
    elseif spec.axes isa Symbol
        String(spec.axes)
    else
        String[String(axis) for axis in spec.axes]
    end
    materialization = isnothing(spec.materialization) ? nothing : Dict{String,Any}(
        "type" => "diameter_indexed",
        "role" => isnothing(spec.materialization.role) ? nothing : String(spec.materialization.role),
        "fill_value" => _encode_value(spec.materialization.fill_value),
    )
    return Dict{String,Any}(
        "name" => String(spec.name),
        "shape" => String(spec.shape),
        "axes" => axes,
        "materialization" => materialization,
    )
end

function _decode_parameter_spec(x, path)
    x = _check_keys(x, ("name", "shape", "axes", "materialization"), path)
    name = _symbol(_required(x, "name", path), "$path.name")
    shape = _symbol(_required(x, "shape", path), "$path.shape")
    shape in (:scalar, :vector, :matrix) || throw(ArgumentError("$path.shape has unsupported shape $(repr(shape))."))
    axes_value = _required(x, "axes", path)
    axes = if axes_value === nothing
        nothing
    elseif axes_value isa AbstractString
        Symbol(axes_value)
    elseif axes_value isa AbstractVector
        length(axes_value) == 2 || throw(ArgumentError("$path.axes must contain exactly two axis names."))
        (_symbol(axes_value[1], "$path.axes[1]"), _symbol(axes_value[2], "$path.axes[2]"))
    else
        throw(ArgumentError("$path.axes must be null, a string, or a two-element array."))
    end
    materialization_value = _required(x, "materialization", path)
    materialization = if materialization_value === nothing
        nothing
    else
        materialization_path = "$path.materialization"
        materialization_value = _check_keys(
            materialization_value, ("type", "role", "fill_value"), materialization_path
        )
        type_name = _string(_required(materialization_value, "type", materialization_path), "$materialization_path.type")
        type_name == "diameter_indexed" || throw(ArgumentError("$materialization_path.type has unsupported materialization $(repr(type_name))."))
        role_value = _required(materialization_value, "role", materialization_path)
        role = isnothing(role_value) ? nothing : _symbol(role_value, "$materialization_path.role")
        fill_value = _decode_value(_required(materialization_value, "fill_value", materialization_path), "$materialization_path.fill_value")
        DiameterIndexedMaterialization(role; fill_value)
    end
    return ParameterSpec(name, shape; axes, materialization)
end

function _encode_default(provider::ConstDefault)
    return Dict{String,Any}("type" => "const", "value" => _encode_value(provider.value))
end
_encode_default(::NoDefault) = Dict{String,Any}("type" => "none")
function _encode_default(provider::FillDefault)
    return Dict{String,Any}("type" => "fill", "value" => _encode_value(provider.value))
end
function _encode_default(provider::DiameterIndexedVectorDefault)
    return Dict{String,Any}(
        "type" => "diameter_indexed_vector",
        "value" => _encode_value(provider.value),
        "role" => String(provider.role),
        "fill_value" => _encode_value(provider.default),
    )
end

function _decode_default(x, path)
    x isa AbstractDict || throw(ArgumentError("$path must be an object."))
    type_name = _string(_required(x, "type", path), "$path.type")
    if type_name == "const"
        _check_keys(x, ("type", "value"), path)
        return ConstDefault(_decode_value(_required(x, "value", path), "$path.value"))
    elseif type_name == "none"
        _check_keys(x, ("type",), path)
        return NoDefault()
    elseif type_name == "fill"
        _check_keys(x, ("type", "value"), path)
        return FillDefault(_decode_value(_required(x, "value", path), "$path.value"))
    elseif type_name == "diameter_indexed_vector"
        _check_keys(x, ("type", "value", "role", "fill_value"), path)
        return DiameterIndexedVectorDefault(
            _decode_value(_required(x, "value", path), "$path.value"),
            _symbol(_required(x, "role", path), "$path.role");
            default=_decode_value(_required(x, "fill_value", path), "$path.fill_value"),
        )
    end
    throw(ArgumentError("$path.type has unsupported default provider $(repr(type_name))."))
end

function _encode_parameter_definitions(definitions)
    return Any[
        Dict{String,Any}(
            "spec" => _encode_parameter_spec(definition.spec),
            "default" => _encode_default(definition.default),
        ) for definition in definitions
    ]
end

function _decode_parameter_definitions(x, path)
    x isa AbstractVector || throw(ArgumentError("$path must be an array."))
    return Tuple(
        begin
            item_path = "$path[$i]"
            item = _check_keys(value, ("spec", "default"), item_path)
            ParameterDefinition(
                _decode_parameter_spec(_required(item, "spec", item_path), "$item_path.spec"),
                _decode_default(_required(item, "default", item_path), "$item_path.default"),
            )
        end for (i, value) in pairs(x)
    )
end

function _encode_interaction_definitions(definitions::NamedTuple)
    return Any[
        Dict{String,Any}(
            "name" => String(name),
            "deriver" => _deriver_id(definition.deriver),
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
        deriver = _decode_deriver(deriver_name, "$item_path.deriver")
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

"""Encode a `ModelRecipe` as a JSON-compatible typed recipe document."""
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
        "scalar_type" => _scalar_type_id(recipe.scalar_type),
    )

    return Dict{String,Any}(
        "schema" => MODEL_RECIPE_SCHEMA,
        "model" => Dict{String,Any}("family" => String(recipe.family)),
        "recipe" => data,
    )
end

"""Decode a typed recipe document into a `ModelRecipe`."""
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
        _decode_scalar_type(recipe["scalar_type"], "Recipe document.recipe.scalar_type"),
    )
end

function decode_recipe(document)
    throw(ArgumentError("Expected an Agate recipe dictionary, got $(typeof(document))."))
end

"""Write a typed recipe document to `path` as pretty-printed JSON."""
function export_recipe(path::AbstractString, recipe::ModelRecipe)
    document = encode_recipe(recipe)
    open(path, "w") do io
        JSON.print(io, document, 4)
        println(io)
    end
    return path
end

"""Read and decode a typed recipe document from `path`."""
import_recipe(path::AbstractString) = decode_recipe(JSON.parsefile(path))
