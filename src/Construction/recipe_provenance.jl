using SHA

const _PROVENANCE_KEYS = ("agate", "provider")
const _PACKAGE_PROVENANCE_KEYS = ("package", "version", "repository", "commit")

function _git_output(cmd::Cmd)
    try
        value = strip(read(pipeline(cmd, stderr=devnull), String))
        return isempty(value) ? nothing : String(value)
    catch
        return nothing
    end
end

function _git_implementation_is_clean(root)
    try
        status = read(
            pipeline(
                `git -C $root status --porcelain --untracked-files=all -- .`,
                stderr=devnull,
            ), String
        )
        return isempty(strip(status))
    catch
        return false
    end
end

_sanitize_repository_url(::Nothing) = nothing
function _sanitize_repository_url(repository::AbstractString)
    url = String(repository)
    occursin("://", url) || return url
    url = replace(url, r"^([A-Za-z][A-Za-z0-9+.-]*://)[^/@]+@" => s"\1")
    return replace(url, r"[?#].*$" => "")
end

function _package_provenance(mod::Module)
    root = Base.moduleroot(mod)
    record = Dict{String,Any}("package" => String(nameof(root)))
    version = Base.pkgversion(root)
    isnothing(version) || (record["version"] = string(version))

    package_root = Base.pkgdir(root)
    isnothing(package_root) && return record
    repository = _sanitize_repository_url(
        _git_output(`git -C $package_root config --get remote.origin.url`)
    )
    commit = _git_implementation_is_clean(package_root) ?
             _git_output(`git -C $package_root rev-parse HEAD`) : nothing
    isnothing(repository) || (record["repository"] = repository)
    isnothing(commit) || (record["commit"] = commit)
    return record
end

function _recipe_provenance(recipe::ModelRecipe)
    agate = Base.moduleroot(@__MODULE__)
    provenance = Dict{String,Any}("agate" => _package_provenance(agate))
    provider = Base.moduleroot(parentmodule(typeof(registered_family(Val(recipe.family)))))
    provider === agate || (provenance["provider"] = _package_provenance(provider))
    return provenance
end

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

function _recipe_hash(family::Symbol, definition_version::VersionNumber, realization)
    content = Dict{String,Any}(
        "schema" => MODEL_RECIPE_SCHEMA,
        "family" => String(family),
        "definition_version" => string(definition_version),
        "realization" => realization,
    )
    return "sha256:" * bytes2hex(sha256(_canonical_json(content)))
end

function _decode_package_provenance(x, path)
    x = _check_keys(x, _PACKAGE_PROVENANCE_KEYS, path)
    _string(_required(x, "package", path), "$path.package")
    for key in ("version", "repository", "commit")
        haskey(x, key) && _string(x[key], "$path.$key")
    end
    return x
end

function _decode_provenance(x, path)
    x = _check_keys(x, _PROVENANCE_KEYS, path)
    _decode_package_provenance(_required(x, "agate", path), "$path.agate")
    haskey(x, "provider") && _decode_package_provenance(x["provider"], "$path.provider")
    return x
end

function _check_recipe_provenance(recipe::ModelRecipe, recorded)
    current = _recipe_provenance(recipe)
    haskey(recorded, "provider") == haskey(current, "provider") ||
        @warn "Model-provider provenance differs from the loaded recipe family."
    for (key, label) in (("agate", "Agate"), ("provider", "Model provider"))
        haskey(recorded, key) && haskey(current, key) || continue
        for field in ("package", "version", "commit")
            haskey(recorded[key], field) && haskey(current[key], field) || continue
            recorded[key][field] == current[key][field] && continue
            @warn "$label $field differs from recipe provenance." recorded=recorded[key][field] loaded=current[key][field]
        end
    end
    return nothing
end
