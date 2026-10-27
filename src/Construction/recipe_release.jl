using Pkg
using TOML

# Keep package-source checks independent of the model-family and Pkg versions.
function _check_release_source(url::AbstractString, name)
    scheme = match(r"^([A-Za-z][A-Za-z0-9+.-]*)://", url)
    protocol = isnothing(scheme) ? nothing : lowercase(scheme[1])
    remote = isnothing(protocol) ? occursin(r"^[^/@:\\]+@[^/:\\]+:.+", url) :
             protocol in ("https", "http", "ssh", "git")
    remote || throw(ArgumentError(
        "Cannot export a portable model release: local source paths for $name ($url)."
    ))
    authority = isnothing(protocol) ? "" : first(split(split(url, "://"; limit=2)[2], r"[/?#]"; limit=2))
    (occursin(r"[?#]", url) ||
     (protocol in ("http", "https") &&
      occursin('@', authority))) && throw(ArgumentError(
        "Cannot export a model release with credentials or query parameters in the source URL for $name."
    ))
    return nothing
end

_release_entries(entries::AbstractVector) = entries
_release_entries(entry) = (entry,)

function _release_project_path()
    project = Base.active_project()
    !isnothing(project) && isfile(project) || throw(ArgumentError(
        "Activate a Julia project with a Project.toml before exporting a model release."
    ))
    return project
end

function _release_manifest_path(project_path)
    # Julia 1.10-1.12 expose this selector; Julia 1.13 also offers Base.active_manifest().
    path = Base.project_file_manifest_path(project_path)
    !isnothing(path) && isfile(path) || throw(ArgumentError(
        "A model release requires a resolved Julia manifest. Run Pkg.instantiate() in the active project."
    ))
    return path
end

function _release_dependencies(manifest)
    format = try
        VersionNumber(string(get(manifest, "manifest_format", "1.0")))
    catch
        nothing
    end
    dependencies = get(manifest, "deps", nothing)
    (!isnothing(format) && format.major == 2 && dependencies isa AbstractDict) || throw(ArgumentError(
        "Model releases require a resolved Julia manifest (format 2)."
    ))
    return dependencies
end

# Verify that the selected manifest describes the modules actually loaded in this session.
function _release_environment(recipe::ModelRecipe, provenance)
    project_path = _release_project_path()
    manifest_path = _release_manifest_path(project_path)
    project, manifest = TOML.parsefile(project_path), TOML.parsefile(manifest_path)
    dependencies = _release_dependencies(manifest)
    haskey(project, "workspace") && throw(ArgumentError(
        "Export model releases from a standalone Julia project rather than a multi-project workspace."
    ))

    # The exported Project must select the copied manifest when instantiated elsewhere.
    selected_name = "Manifest.toml"
    if haskey(project, "manifest")
        explicit = project["manifest"]
        explicit isa AbstractString || throw(ArgumentError("The project has an invalid manifest path."))
        selected_name = normpath(explicit)
        (!isabspath(selected_name) && selected_name != "." &&
         !(".." in splitpath(selected_name))) ||
            throw(ArgumentError("The project specifies a nonportable manifest path."))
    end

    local_packages = String[]
    for (name, source) in get(project, "sources", Dict())
        source isa AbstractDict || continue
        haskey(source, "path") && push!(local_packages, name)
        haskey(source, "url") && _check_release_source(source["url"], name)
    end
    for (name, registry) in get(manifest, "registries", Dict())
        haskey(registry, "url") && _check_release_source(registry["url"], name)
    end
    for (name, entries) in dependencies, entry in _release_entries(entries)
        haskey(entry, "path") && push!(local_packages, name)
        if haskey(entry, "repo-url")
            _check_release_source(entry["repo-url"], name)
            occursin(r"^[0-9a-fA-F]{40}$", get(entry, "git-tree-sha1", "")) || throw(ArgumentError(
                "The release manifest does not pin the source tree for $name."
            ))
        end
    end
    isempty(local_packages) || throw(ArgumentError(
        "Cannot export a portable model release: local source paths for " *
        join(sort!(unique!(local_packages)), ", ") * "."
    ))

    project_deps = get(project, "deps", Dict())
    for (name, uuid) in project_deps
        count(entry -> get(entry, "uuid", nothing) == uuid,
              _release_entries(get(dependencies, name, ()))) == 1 || throw(ArgumentError(
            "The release manifest does not resolve direct dependency $name ($uuid)."
        ))
    end
    # Pkg resolves version, source tree, and source path from the selected manifest.
    resolved = Pkg.dependencies()
    tracked_paths = sort!(String[info.name for info in values(resolved) if info.is_tracking_path])
    isempty(tracked_paths) || throw(ArgumentError(
        "Cannot export a portable model release: local source paths for " * join(tracked_paths, ", ") * "."
    ))
    for (label, mod) in (("agate", Base.moduleroot(@__MODULE__)),
                         ("provider", Base.moduleroot(parentmodule(typeof(registered_family(Val(recipe.family)))))))
        label == "provider" && mod === Base.moduleroot(@__MODULE__) && continue
        package_id = Base.PkgId(mod)
        name = String(package_id.name)
        get(project_deps, name, nothing) == string(package_id.uuid) || throw(ArgumentError(
            "The release environment must include the loaded $name package as a direct dependency."
        ))
        info = get(resolved, package_id.uuid, nothing)
        isnothing(info) && throw(ArgumentError("The release environment cannot resolve loaded package $name."))
        info.version == Base.pkgversion(mod) || throw(ArgumentError(
            "The release manifest does not match the loaded $name version."
        ))
        !isnothing(info.tree_hash) || throw(ArgumentError(
            "The release manifest does not pin the loaded $name source tree."
        ))
        loaded_source = Base.pkgdir(mod)
        !isnothing(loaded_source) && isdir(info.source) && isdir(loaded_source) &&
            realpath(info.source) == realpath(loaded_source) || throw(ArgumentError(
            "The loaded $name source differs from the release environment. Restart Julia using the release project."
        ))
        record = get(provenance, label, nothing)
        !isnothing(record) && record["package"] == name || throw(ArgumentError(
            "The loaded package identity for $name differs from recipe provenance."
        ))
    end
    for name in ("LocalPreferences.toml", "JuliaLocalPreferences.toml")
        isfile(joinpath(dirname(project_path), name)) &&
            @warn "Local Julia preferences are not included in the model release; document any experiment-specific settings." file=name
    end
    return (project_path, manifest_path, selected_name)
end

"""Export a recipe and its resolved, portable Julia environment to a release directory.

The environment must be the active Julia project, with installed package sources matching
those already loaded. The release contains the selected manifest, checked recipe, and project.
"""
function export_model_release(directory::AbstractString, recipe::ModelRecipe)
    document = encode_recipe(recipe)
    decode_recipe(document) # Validate the same scientific contract as import.
    project, manifest, manifest_name = _release_environment(recipe, document["provenance"])
    destination = abspath(directory)
    if ispath(destination)
        isdir(destination) && !islink(destination) && isempty(readdir(destination)) ||
            throw(ArgumentError("Model release destination must be a new or empty directory: $directory"))
    end
    mkpath(dirname(destination))
    mktempdir(dirname(destination); prefix=".agate-release-") do staging
        open(joinpath(staging, "recipe.json"), "w") do io
            JSON.json(io, document; pretty=4, omit_empty=false)
            println(io)
        end
        cp(project, joinpath(staging, "Project.toml"); follow_symlinks=true)
        target_manifest = joinpath(staging, manifest_name)
        mkpath(dirname(target_manifest))
        cp(manifest, target_manifest; follow_symlinks=true)
        isdir(destination) && isempty(readdir(destination)) && rm(destination)
        mv(staging, destination)
    end
    return directory
end
