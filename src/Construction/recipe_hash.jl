using SHA

function _canonical_json(x)
    if x isa AbstractDict
        keys_sorted = sort!(String[String(k) for k in keys(x)])
        items = (JSON.json(k) * ":" * _canonical_json(x[k]) for k in keys_sorted)
        return "{" * join(items, ",") * "}"
    elseif x isa AbstractVector
        return "[" * join((_canonical_json(v) for v in x), ",") * "]"
    end
    return JSON.json(x)
end

function _recipe_hash(family::Symbol, realization)
    identity = Dict{String,Any}(
        "schema" => MODEL_RECIPE_SCHEMA,
        "family" => String(family),
        "realization" => realization,
    )
    return "sha256:" * bytes2hex(sha256(_canonical_json(identity)))
end
