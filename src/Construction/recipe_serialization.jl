using JSON

using ..Configuration:
    PFTSpecification, DiameterListSpecification, DiameterRangeSpecification, Population, Pool,
    currency, states, state_currency, size_structure, normalize_diameters

using ..ModelFamilies: default_components, default_processes
using ..Parameters: parameter_definitions

using ..Processes:
    AbstractProcess, AbstractFormulation, AbstractFactor, FactorDriver, FactorComponent,
    ModelDefinition, normalize_model, parameter_bindings, process_id, process_kind, factor_kind, formulation,
    formulation_tag, formulation_recipe_fields, factors, factor_inputs, factor_children,
    factor_child_path, participants, drivers, rate_axes,
    process_routing, process_products, product_path, product_recipe_key, process_stoichiometry,
    Growth, Nutrients, Products, ProductRouting,
    FixedStoichiometry

using ..Library.Allometry:
    ConstantParam,
    AllometricParam,
    allometric_relationship_identifier,
    allometric_relationship_from_identifier

const PROCESS_MODEL_RECIPE_SCHEMA = "agate.model_recipe.v0.7"
const _RECIPE_DOCUMENT_KEYS = ("schema", "model", "provenance", "recipe", "recipe_hash")
const _RECIPE_MODEL_KEYS = ("family",)
const _SUPPORTED_SPLITTING = (:linear_splitting, :log_splitting)

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

function _symbol(x, path)
    return Symbol(_string(x, path))
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

_json_value(x::Nothing) = nothing
_json_value(x::Bool) = x
_json_value(x::Int) = x
_json_value(x::AbstractFloat) = _finite_float(x)
_json_value(x::AbstractString) = String(x)

function _json_value(x::AbstractDict)
    all(key -> key isa AbstractString, keys(x)) || throw(
        ArgumentError("Recipe JSON objects must use string keys.")
    )
    return Dict{String,Any}(String(key) => _json_value(value) for (key, value) in pairs(x))
end

function _json_value(x::AbstractVector)
    return Any[_json_value(value) for value in x]
end

function _json_value(x)
    throw(ArgumentError("Recipe JSON data has unsupported value of type $(typeof(x))."))
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

_size_structure_data(spec::DiameterListSpecification) = Dict{String,Any}(
    "diameters" => [_finite_float(float(value)) for value in spec.diameters]
)
_size_structure_data(spec::DiameterRangeSpecification) = Dict{String,Any}(
    "n" => spec.n,
    "min_esd" => _finite_float(float(spec.min_diameter)),
    "max_esd" => _finite_float(float(spec.max_diameter)),
    "splitting" => String(spec.splitting),
)
_size_structure_data(spec) = _size_structure_data(normalize_diameters(spec).specification)

function _component_recipe_data(name::Symbol, component::Population, recipe::ProcessModelRecipe)
    groups = getproperty(recipe.population_groups, name)
    data = Dict{String,Any}(
        "kind" => "population",
        "states" => Dict(
            String(state) => String(state_currency(component, state)) for state in keys(states(component))
        ),
        "size_structure" => nothing,
    )
    if length(groups) == 1 && only(groups) === name
        data["size_structure"] = _size_structure_data(getproperty(recipe.community, name).diameters)
    else
        data["subgroups"] = Dict{String,Any}(
            String(group) => _size_structure_data(getproperty(recipe.community, group).diameters)
            for group in groups
        )
    end
    return data
end

function _component_recipe_data(::Symbol, component::Pool, ::ProcessModelRecipe)
    structure = size_structure(component)
    return Dict{String,Any}(
        "kind" => "pool",
        "currency" => String(currency(component)),
        "size_structure" => isnothing(structure) ? nothing : _size_structure_data(structure),
    )
end

_recipe_science_value(x::Nothing) = nothing
_recipe_science_value(x::Bool) = x
_recipe_science_value(x::Integer) = x
_recipe_science_value(x::AbstractFloat) = _finite_float(x)
_recipe_science_value(x::AbstractString) = String(x)
_recipe_science_value(x::Symbol) = String(x)
_recipe_science_value(x::Tuple) = Any[_recipe_science_value(value) for value in x]
_recipe_science_value(x::NamedTuple) = Dict{String,Any}(
    String(name) => _recipe_science_value(value) for (name, value) in pairs(x)
)

function _participants_recipe_data(named)
    return Dict{String,Any}(
        String(role) => _recipe_science_value(length(values) == 1 ? only(values) : values)
        for (role, values) in pairs(participants(named))
    )
end

function _validated_formulation_recipe_fields(formulation_value::AbstractFormulation)
    fields = formulation_recipe_fields(formulation_value)
    fields isa NamedTuple || throw(ArgumentError(
        "formulation_recipe_fields for $(typeof(formulation_value)) must return a NamedTuple",
    ))
    return fields
end

function _attach_formulation_recipe_fields!(data::Dict{String,Any}, formulation_value)
    formulation_value isa AbstractFormulation || return data
    fields = _validated_formulation_recipe_fields(formulation_value)
    isempty(fields) || (data["formulation_fields"] = _recipe_science_value(fields))
    return data
end

function _node_bindings_recipe_data(normalized, process::Symbol, path::Tuple)
    matches = (
        binding for binding in parameter_bindings(normalized)
        if binding.process === process && binding.path == path
    )
    result = Dict{String,Any}()
    for binding in matches
        slot = String(binding.slot)
        if isempty(binding.qualifier)
            haskey(result, slot) && throw(ArgumentError(
                "recipe serialization found repeated unqualified binding for process :$process path $path slot :$(binding.slot)",
            ))
            result[slot] = String(binding.parameter)
            continue
        end

        length(binding.qualifier) == 1 || throw(ArgumentError(
            "recipe serialization supports one qualifier dimension per parameter slot",
        ))
        qualified = get!(result, slot, Dict{String,Any}())
        qualified isa Dict || throw(ArgumentError(
            "recipe serialization found mixed qualified and unqualified bindings for process :$process path $path slot :$(binding.slot)",
        ))
        qualified[String(only(values(binding.qualifier)))] = String(binding.parameter)
    end
    return isempty(result) ? nothing : result
end

function _attach_node_bindings!(data::Dict{String,Any}, normalized, process::Symbol, path::Tuple)
    bindings = _node_bindings_recipe_data(normalized, process, path)
    isnothing(bindings) || (data["bindings"] = bindings)
    return data
end

_factor_input_fields(input::FactorDriver) = (drivers=(driver=input.identity,),)
_factor_input_fields(input::FactorComponent) = (participants=(resource=input.component,),)
_factor_children_key(::AbstractFactor) = :factors
_factor_children_key(::Nutrients) = :responses

function _factor_recipe_data(
    factor::AbstractFactor, normalized, process::Symbol, path::Tuple
)
    factor_formulation = formulation(factor)
    data = Dict{String,Any}(
        "kind" => String(factor_kind(factor)),
        "formulation" => String(formulation_tag(factor_formulation)),
    )
    _attach_formulation_recipe_fields!(data, factor_formulation)

    for input in factor_inputs(factor), (name, value) in pairs(_factor_input_fields(input))
        data[String(name)] = _recipe_science_value(value)
    end
    _attach_node_bindings!(data, normalized, process, path)

    children = factor_children(factor)
    if !isempty(children)
        child_data = Dict{String,Any}()
        for (name, child) in pairs(children)
            child_path = factor_child_path(path, factor, name)
            child_data[String(name)] = _factor_recipe_data(child, normalized, process, child_path)
        end
        data[String(_factor_children_key(factor))] = child_data
    end
    return data
end

function _stoichiometry_recipe_data(
    stoichiometry::FixedStoichiometry, normalized, process::Symbol, path::Tuple
)
    data = Dict{String,Any}(
        "kind" => String(formulation_tag(stoichiometry)),
        "reference" => String(stoichiometry.reference),
    )
    _attach_node_bindings!(data, normalized, process, path)
    return data
end

function _products_recipe_data(
    products::Products, _normalized, _process::Symbol, _path::Tuple
)
    if length(products.targets) == 1
        return String(only(values(products.targets)))
    end
    return Dict{String,Any}(
        "targets" => _recipe_science_value(products.targets),
        "fractions" => _recipe_science_value(products.fractions),
    )
end

function _routing_recipe_data(
    routing::ProductRouting, normalized, process::Symbol
)
    data = Dict{String,Any}(
        "kind" => "product_routing",
        "formulation" => String(formulation_tag(formulation(routing))),
        "pools" => _recipe_science_value(routing.pools),
        "stoichiometry" => _stoichiometry_recipe_data(
            routing.stoichiometry, normalized, process, (:routing, :stoichiometry)
        ),
    )
    _attach_node_bindings!(data, normalized, process, (:routing,))
    return data
end

_process_formulation_tag(::Growth) = nothing
_process_formulation_tag(process::AbstractProcess) = formulation_tag(formulation(process))

function _process_recipe_data(named, normalized)
    process = named.process
    process_name = process_id(named)
    data = Dict{String,Any}(
        "kind" => String(process_kind(named)),
        "participants" => _participants_recipe_data(named),
        "rate_axes" => _recipe_science_value(rate_axes(named)),
    )

    formulation_name = _process_formulation_tag(process)
    if !isnothing(formulation_name)
        process_formulation = formulation(process)
        data["formulation"] = String(formulation_name)
        _attach_formulation_recipe_fields!(data, process_formulation)
    end
    _attach_node_bindings!(data, normalized, process_name, ())

    process_factors = factors(named)
    if !isempty(process_factors)
        factor_data = Dict{String,Any}()
        for (name, factor) in pairs(process_factors)
            factor_data[String(name)] = _factor_recipe_data(
                factor, normalized, process_name, (:factors, name)
            )
        end
        data["factors"] = factor_data
    else
        process_drivers = drivers(named)
        isempty(process_drivers) || (data["drivers"] = _recipe_science_value(process_drivers))
    end

    stoichiometry = process_stoichiometry(process)
    isnothing(stoichiometry) || (
        data["stoichiometry"] = _stoichiometry_recipe_data(
            stoichiometry, normalized, process_name, (:stoichiometry,)
        )
    )
    products = process_products(process)
    if !isnothing(products)
        path = product_path(process)
        data[String(product_recipe_key(process))] = _products_recipe_data(
            products, normalized, process_name, path
        )
    else
        routing = process_routing(process)
        isnothing(routing) || (
            data["routing"] = _routing_recipe_data(routing, normalized, process_name)
        )
    end
    return data
end

function _encode_process_recipe_data(recipe::ProcessModelRecipe)
    normalized = normalize_model(ModelDefinition(;
        components=recipe.components,
        processes=recipe.processes,
        parameters=parameter_definitions(registered_family(Val(recipe.family))),
    ))
    components = Dict{String,Any}(
        String(name) => _component_recipe_data(name, getproperty(recipe.components, name), recipe)
        for name in keys(recipe.components)
    )
    processes = Dict{String,Any}(
        String(name) => _process_recipe_data(getproperty(normalized.processes, name), normalized)
        for name in keys(normalized.processes)
    )
    realization = Dict{String,Any}(
        "community" => _encode_community(recipe.community),
        "population_groups" => _encode_value(recipe.population_groups),
        "parameter_overrides" => _encode_value(recipe.parameter_overrides),
        "sinking_tracers" => isnothing(recipe.sinking_tracers) ? nothing : _encode_value(recipe.sinking_tracers),
        "open_bottom" => recipe.open_bottom,
    )
    return Dict{String,Any}(
        "components" => components,
        "processes" => processes,
        "realization" => realization,
    )
end

"""Encode a recipe with its scientific hash and available package provenance."""
function encode_recipe(recipe::ProcessModelRecipe)
    data = _encode_process_recipe_data(recipe)
    return _json_value(Dict{String,Any}(
        "schema" => PROCESS_MODEL_RECIPE_SCHEMA,
        "model" => Dict{String,Any}("family" => String(recipe.family)),
        "provenance" => _recipe_provenance(recipe),
        "recipe" => data,
        "recipe_hash" => _recipe_hash(recipe.family, data),
    ))
end

"""Decode a recipe document, verifying its hash and checking package provenance."""
const _PROCESS_RECIPE_KEYS = ("components", "processes", "realization")
const _PROCESS_REALIZATION_KEYS = (
    "community",
    "population_groups",
    "parameter_overrides",
    "sinking_tracers",
    "open_bottom",
)

function _decode_process_model_recipe(document::AbstractDict)
    model = _check_keys(
        _required(document, "model", "Recipe document"), _RECIPE_MODEL_KEYS, "Recipe document.model"
    )
    family = _symbol(_required(model, "family", "Recipe document.model"), "Recipe document.model.family")
    provenance = _decode_provenance(
        _required(document, "provenance", "Recipe document"), "Recipe document.provenance"
    )
    recorded_hash = _string(
        _required(document, "recipe_hash", "Recipe document"), "Recipe document.recipe_hash"
    )
    recipe_data = _complete_object(
        _required(document, "recipe", "Recipe document"), _PROCESS_RECIPE_KEYS, "Recipe document.recipe"
    )
    recorded_hash == _recipe_hash(family, recipe_data) || throw(
        ArgumentError("Recipe document.recipe_hash does not match the serialized recipe content.")
    )
    model_family = registered_family(Val(family))
    realization = _complete_object(
        recipe_data["realization"], _PROCESS_REALIZATION_KEYS, "Recipe document.recipe.realization"
    )
    community = _decode_community(realization["community"], "Recipe document.recipe.realization.community")
    population_groups = _decode_value(
        realization["population_groups"], "Recipe document.recipe.realization.population_groups"
    )
    population_groups isa NamedTuple || throw(
        ArgumentError("Recipe document.recipe.realization.population_groups must decode to a NamedTuple.")
    )
    all(groups -> groups isa Tuple && all(group -> group isa Symbol, groups), values(population_groups)) || throw(
        ArgumentError("Recipe population_groups values must be tuples of Symbols.")
    )
    parameter_overrides = _decode_value(
        realization["parameter_overrides"], "Recipe document.recipe.realization.parameter_overrides"
    )
    parameter_overrides isa NamedTuple || throw(
        ArgumentError("Recipe parameter_overrides must decode to a NamedTuple.")
    )
    sinking_tracers = isnothing(realization["sinking_tracers"]) ? nothing : _decode_value(
        realization["sinking_tracers"], "Recipe document.recipe.realization.sinking_tracers"
    )
    !isnothing(sinking_tracers) && !(sinking_tracers isa NamedTuple) && throw(
        ArgumentError("Recipe sinking_tracers must decode to a NamedTuple or null.")
    )

    decoded = ProcessModelRecipe(
        family,
        deepcopy(default_components(model_family)),
        deepcopy(default_processes(model_family)),
        population_groups,
        community,
        parameter_overrides,
        sinking_tracers,
        _boolean(realization["open_bottom"], "Recipe document.recipe.realization.open_bottom"),
    )

    expected_science = _encode_process_recipe_data(decoded)
    recorded_hash == _recipe_hash(family, expected_science) || throw(
        ArgumentError("Recipe document does not match the loaded model family contract.")
    )
    _check_recipe_provenance(decoded, provenance)
    return decoded
end

"""Decode and validate an Agate recipe document for the current recipe schema."""
function decode_recipe(document::AbstractDict)
    document = _check_keys(document, _RECIPE_DOCUMENT_KEYS, "Recipe document")
    schema = _string(_required(document, "schema", "Recipe document"), "Recipe document.schema")
    schema == PROCESS_MODEL_RECIPE_SCHEMA || throw(
        ArgumentError(
            "Unsupported Agate recipe schema $(repr(schema)); supported schema is " *
            "$(repr(PROCESS_MODEL_RECIPE_SCHEMA))."
        )
    )
    return _decode_process_model_recipe(document)
end
function decode_recipe(document)
    throw(ArgumentError("Expected an Agate recipe dictionary, got $(typeof(document))."))
end

"""Write a recipe document to `path` as pretty-printed JSON."""
function export_recipe(path::AbstractString, recipe::ProcessModelRecipe)
    document = encode_recipe(recipe)
    open(path, "w") do io
        JSON.json(io, document; pretty=4, omit_empty=false)
        println(io)
    end
    return path
end

"""Read and decode a recipe document from `path`."""
import_recipe(path::AbstractString) = decode_recipe(JSON.parsefile(path))
