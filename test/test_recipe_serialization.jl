using Agate.Construction: decode_recipe, encode_recipe, export_recipe, import_recipe
using Agate.Models: NiPiZD
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
    _, default_recipe = NiPiZD.construct_with_recipe()
    @test decode_recipe(encode_recipe(default_recipe)) == default_recipe

    _, recipe, manifest = NiPiZD.construct_with_manifest(;
        authored_nipizd_inputs(Float32)...
    )

    encoded = encode_recipe(recipe)
    decoded = decode_recipe(encoded)
    _, decoded_manifest = NiPiZD.construct_with_manifest(decoded)

    @test encoded["schema"] == "agate.model_recipe.v1"
    @test Set(keys(encoded)) == Set(("schema", "model", "provenance", "recipe", "recipe_hash"))
    @test !haskey(encoded["recipe"], "parameter_definitions")
    @test !haskey(encoded["recipe"], "interaction_definitions")
    @test encoded["provenance"]["agate"]["package"] == "Agate"
    @test encoded["provenance"]["agate"]["version"] == "0.10.0"
    @test startswith(encoded["recipe_hash"], "sha256:")
    @test length(encoded["recipe_hash"]) == 71
    @test encode_recipe(recipe)["recipe_hash"] == encoded["recipe_hash"]
    @test decoded == recipe
    @test decoded_manifest == manifest
    @test decoded_manifest.sinking_tracers.D isa Float32

    mktemp() do path, io
        close(io)
        @test export_recipe(path, recipe) == path
        @test import_recipe(path) == recipe
    end

    _, nonfinite_recipe = NiPiZD.construct_with_recipe(;
        scalar_type=Float32, palatability_matrix=Float32[NaN 0.9; 0.3 0.7]
    )
    @test_throws ArgumentError encode_recipe(nonfinite_recipe)

    version_mismatch = modified(
        encoded, x -> (x["provenance"]["agate"]["version"] = "0.9.0")
    )
    warned = @test_logs (:warn, r"Agate version differs") decode_recipe(version_mismatch)
    @test warned == recipe

    invalid_documents = (
        modified(encoded, x -> (x["extra"] = true)),
        modified(encoded, x -> delete!(x["recipe"], "open_bottom")),
        modified(encoded, x -> (x["schema"] = "agate.model_recipe.v2")),
        modified(encoded, x -> (x["model"]["family"] = "UnknownModel")),
        modified(encoded, x -> (x["recipe"]["scalar_type"] = "Float16")),
        modified(encoded, x -> (x["recipe"]["open_bottom"] = !x["recipe"]["open_bottom"])),
        modified(encoded, x -> (
            encoded_named_value(x["recipe"]["parameter_overrides"], :linear_mortality)["law"] = "unknown"
        )),
        modified(encoded, x -> (
            encoded_named_value(x["recipe"]["interaction_overrides"], :palatability_matrix)[2] = [0.3]
        )),
    )

    for document in invalid_documents
        @test_throws ArgumentError decode_recipe(document)
    end
    @test_throws ArgumentError decode_recipe("not a recipe document")
end
