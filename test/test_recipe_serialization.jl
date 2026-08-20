using Agate.Construction: decode_recipe, encode_recipe, export_recipe, import_recipe
using Agate.Models: NiPiZD
using JSON
using OceanBioME: BoxModelGrid
using Test

function encoded_named_value(encoded, name)
    entry = only(item for item in encoded["entries"] if item["name"] == String(name))
    return entry["value"]
end

function modified(document, mutation)
    copy = deepcopy(document)
    mutation(copy)
    return copy
end

@testset "NiPiZD recipe serialization" begin
    _, default_recipe = NiPiZD.construct_plus_recipe()
    default_encoded = encode_recipe(default_recipe)
    @test decode_recipe(default_encoded) == default_recipe
    @test default_encoded["schema"] == "agate.model_recipe.v3"

    target = JSON.parsefile(
        joinpath(@__DIR__, "reference", "nipizd_model_recipe_v3_target.json")
    )["recipe"]
    for key in ("components", "processes", "parameter_bindings")
        @test default_encoded["recipe"][key] == target[key]
    end

    bgc, recipe = NiPiZD.construct_plus_recipe(; authored_nipizd_inputs(Float32)...)
    manifest = nipizd_manifest(recipe; scalar_type=Float32)

    encoded = encode_recipe(recipe)
    realization = encoded["recipe"]["realization"]
    recipe_hash = encoded["recipe_hash"]
    bgc.parameters.palatability_matrix[1, 1] = 0f0
    @test recipe.interaction_overrides.palatability_matrix[1, 1] == 0.8f0
    @test encode_recipe(recipe)["recipe_hash"] == recipe_hash
    decoded = decode_recipe(encoded)
    decoded_manifest = nipizd_manifest(decoded; scalar_type=Float32)

    @test Set(keys(encoded)) == Set(("schema", "model", "provenance", "recipe", "recipe_hash"))
    @test Set(keys(encoded["recipe"])) ==
        Set(("components", "processes", "parameter_bindings", "realization"))
    @test !haskey(realization, "ecological_roles")
    @test !haskey(realization, "interaction_roles")
    @test !haskey(realization, "parameter_roles")
    @test !haskey(realization, "auxiliary_fields")
    @test !haskey(realization, "scalar_type")
    @test encoded["provenance"]["agate"]["package"] == "Agate"
    @test encoded["provenance"]["agate"]["version"] == string(Base.pkgversion(Agate))
    @test startswith(recipe_hash, "sha256:")
    @test length(recipe_hash) == 71
    @test decoded == recipe
    @test decoded_manifest == manifest
    @test decoded_manifest.sinking_tracers.D isa Float32

    sinking = (D=2.5 / 86400,)
    _, recipe_f32 = NiPiZD.construct_plus_recipe(;
        grid=BoxModelGrid(Float32), sinking_tracers=sinking, open_bottom=false
    )
    _, recipe_f64 = NiPiZD.construct_plus_recipe(;
        grid=BoxModelGrid(Float64), sinking_tracers=sinking, open_bottom=false
    )
    @test recipe_f32 == recipe_f64
    @test encode_recipe(recipe_f32)["recipe_hash"] == encode_recipe(recipe_f64)["recipe_hash"]

    manifest_f32 = nipizd_manifest(recipe_f32; grid=BoxModelGrid(Float32))
    manifest_f64 = nipizd_manifest(recipe_f32; grid=BoxModelGrid(Float64))
    @test manifest_f32.scalar_type === Float32
    @test manifest_f64.scalar_type === Float64
    @test manifest_f32.sinking_tracers.D isa Float32
    @test manifest_f64.sinking_tracers.D isa Float64

    mktemp() do path, io
        close(io)
        @test export_recipe(path, recipe) == path
        exported = JSON.parsefile(path)
        bindings = exported["recipe"]["parameter_bindings"]
        empty_arrays = (
            bindings["alpha"]["provides"][1]["path"],
            bindings["mortality_export_fraction"]["provides"][1]["axes"],
        )
        @test all(value -> value isa AbstractVector && isempty(value), empty_arrays)
        @test import_recipe(path) == recipe
    end

    empty_object_path = modified(encoded, x ->
        (x["recipe"]["parameter_bindings"]["alpha"]["provides"][1]["path"] = Dict{String,Any}())
    )
    @test decode_recipe(empty_object_path) == recipe

    _, nonfinite_recipe = NiPiZD.construct_plus_recipe(;
        scalar_type=Float32, palatability_matrix=Float32[NaN 0.9; 0.3 0.7]
    )
    @test_throws ArgumentError encode_recipe(nonfinite_recipe)

    version_mismatch = modified(
        encoded, x -> (x["provenance"]["agate"]["version"] = "0.9.0")
    )
    warned = @test_logs (:warn, r"Agate version differs") decode_recipe(version_mismatch)
    @test warned == recipe

    v1_document = modified(encoded, x -> (x["schema"] = "agate.model_recipe.v1"))
    v1_error = try
        decode_recipe(v1_document)
        nothing
    catch err
        err
    end
    @test v1_error isa ArgumentError
    message = sprint(showerror, v1_error)
    @test occursin("Unsupported Agate recipe schema", message)
    @test occursin("agate.model_recipe.v2", message)
    @test occursin("agate.model_recipe.v3", message)

    invalid_documents = (
        modified(encoded, x -> (x["extra"] = true)),
        modified(encoded, x -> delete!(x["recipe"]["realization"], "open_bottom")),
        modified(encoded, x -> (x["model"]["family"] = "UnknownModel")),
        modified(encoded, x -> (x["recipe"]["realization"]["scalar_type"] = "Float32")),
        modified(encoded, x -> (x["recipe"]["realization"]["open_bottom"] =
            !x["recipe"]["realization"]["open_bottom"])),
        modified(encoded, x -> (x["recipe"]["processes"]["growth_P"]["formulation"] = "unknown")),
        modified(encoded, x -> (
            encoded_named_value(
                x["recipe"]["realization"]["parameter_overrides"], :linear_mortality
            )["law"] = "unknown"
        )),
        modified(encoded, x -> (
            encoded_named_value(
                x["recipe"]["realization"]["interaction_overrides"], :palatability_matrix
            )[2] = [0.3]
        )),
    )

    for document in invalid_documents
        @test_throws ArgumentError decode_recipe(document)
    end
    @test_throws ArgumentError decode_recipe("not a recipe document")
end

@testset "Legacy v2 recipe round-trip" begin
    factory = Agate.Models.NiPiZD.NiPiZDFactory()
    recipe = Agate.Construction.capture_model_recipe(
        factory;
        community=Agate.Factories.default_community(factory),
        ecological_roles=(phytoplankton=(:P,), zooplankton=(:Z,)),
        interaction_roles=(consumers=(:Z,), prey=(:P,)),
        parameter_roles=(producers=(:P,), consumers=(:Z,)),
        auxiliary_fields=(:PAR,),
    )
    encoded = encode_recipe(recipe)

    @test encoded["schema"] == "agate.model_recipe.v2"
    @test decode_recipe(encoded) == recipe
end
