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

function finalize_construction_manifest(context)
    model = context.model
    resolved = context.resolved

    return Dict{String,Any}(
        "schema" => "agate.construction_manifest.v1",
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
            "scalar_type" => resolved.scalar_type,
            "architecture" => resolved.architecture,
            "has_sinking_velocities" => resolved.has_sinking_velocities,
        ),
    )
end
