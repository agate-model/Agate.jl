using TOML

# Verify that a release environment refers to portable, pinned package sources.
function _release_environment(project_path, recipe::ModelRecipe)
    project_path isa AbstractString && isfile(project_path) || throw(ArgumentError(
        "A model release requires an existing Julia Project.toml. Activate a project with Agate installed."
    ))
    manifest_path = joinpath(dirname(project_path), "Manifest.toml")
    isfile(manifest_path) || throw(ArgumentError(
        "A model release requires a Manifest.toml. Run Pkg.instantiate() in the active project."
    ))
    project, manifest = TOML.parsefile(project_path), TOML.parsefile(manifest_path)
    dependencies = get(manifest, "deps", nothing)
    (get(manifest, "manifest_format", nothing) == "2.0" && dependencies isa AbstractDict) || throw(ArgumentError(
        "Model releases require a resolved Julia Manifest.toml (format 2)."
    ))

    local_packages = String[]
    unpinned_packages = String[]
    for (name, source) in get(project, "sources", Dict())
        source isa AbstractDict || continue
        haskey(source, "path") && push!(local_packages, name)
        haskey(source, "url") && !occursin(r"^[0-9a-fA-F]{40}$", get(source, "rev", "")) &&
            push!(unpinned_packages, name)
    end
    for (name, entries) in dependencies
        for entry in (entries isa AbstractVector ? entries : (entries,))
            haskey(entry, "path") && push!(local_packages, name)
            if haskey(entry, "repo-url")
                url = entry["repo-url"]
                (startswith(url, "file:") || startswith(url, "/")) &&
                    push!(local_packages, name)
                (!occursin(r"^[0-9a-fA-F]{40}$", get(entry, "repo-rev", "")) ||
                 !haskey(entry, "git-tree-sha1")) && push!(unpinned_packages, name)
            end
        end
    end
    isempty(local_packages) || throw(ArgumentError(
        "Cannot export a portable model release: local source paths for " *
        join(sort!(unique!(local_packages)), ", ") *
        ". Install published package versions or commit-pinned Git dependencies in a separate Julia project."
    ))
    isempty(unpinned_packages) || throw(ArgumentError(
        "Cannot export a portable model release: Git dependencies without pinned commits: " *
        join(sort!(unique!(unpinned_packages)), ", ") * "."
    ))

    project_deps = get(project, "deps", Dict())
    for record in values(_recipe_provenance(recipe))
        name = record["package"]
        haskey(project_deps, name) && haskey(dependencies, name) || throw(ArgumentError(
            "The release environment must include $name as a direct dependency."
        ))
        entries = dependencies[name]
        matches = [entry for entry in (entries isa AbstractVector ? entries : (entries,))
                   if get(entry, "uuid", nothing) == project_deps[name]]
        length(matches) == 1 || throw(ArgumentError(
            "The release manifest does not resolve the requested $name dependency."
        ))
        haskey(only(matches), "git-tree-sha1") || throw(ArgumentError(
            "The release manifest does not pin the source tree for $name."
        ))
        expected = get(record, "version", nothing)
        (isnothing(expected) || get(only(matches), "version", nothing) == expected) ||
            throw(ArgumentError("The release manifest does not match the loaded $name version."))
    end
    return manifest_path
end

"""Export a recipe and a checked Julia environment to a release directory.

The active Julia project must have a resolved manifest containing Agate and the
model provider, with portable dependencies (registered packages or Git sources
pinned to commits). The directory must be new or empty. This validates the
package references; actual availability is checked when the release is instantiated.
"""
function export_model_release(directory::AbstractString, recipe::ModelRecipe;
                              project_path=Base.active_project())
    manifest_path = _release_environment(project_path, recipe)
    if ispath(directory)
        isdir(directory) && isempty(readdir(directory)) || throw(ArgumentError(
            "Model release destination must be a new or empty directory: $directory"
        ))
    end
    # Validate the recipe before creating release files.
    encode_recipe(recipe)
    mkpath(directory)
    export_recipe(joinpath(directory, "recipe.json"), recipe)
    cp(project_path, joinpath(directory, "Project.toml"))
    cp(manifest_path, joinpath(directory, "Manifest.toml"))
    return directory
end
