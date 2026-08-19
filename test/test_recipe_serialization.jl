using Agate.Construction: decode_recipe, encode_recipe, export_recipe, import_recipe
using Agate.Models: NiPiZD
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
    @test decode_recipe(encode_recipe(default_recipe)) == default_recipe

    _, recipe = NiPiZD.construct_plus_recipe(; authored_nipizd_inputs(Float32)...)
    manifest = nipizd_manifest(recipe; scalar_type=Float32)

    encoded = encode_recipe(recipe)
    decoded = decode_recipe(encoded)
    decoded_manifest = nipizd_manifest(decoded; scalar_type=Float32)

    @test encoded["schema"] == "agate.model_recipe.v2"
    @test Set(keys(encoded)) == Set(("schema", "model", "provenance", "recipe", "recipe_hash"))
    @test !haskey(encoded["recipe"], "parameter_definitions")
    @test !haskey(encoded["recipe"], "interaction_definitions")
    @test !haskey(encoded["recipe"], "scalar_type")
    @test encoded["provenance"]["agate"]["package"] == "Agate"
    @test encoded["provenance"]["agate"]["version"] == string(Base.pkgversion(Agate))
    @test startswith(encoded["recipe_hash"], "sha256:")
    @test length(encoded["recipe_hash"]) == 71
    @test encode_recipe(recipe)["recipe_hash"] == encoded["recipe_hash"]
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
    encoded_f32 = encode_recipe(recipe_f32)
    encoded_f64 = encode_recipe(recipe_f64)
    @test recipe_f32 == recipe_f64
    @test encoded_f32["recipe_hash"] == encoded_f64["recipe_hash"]

    manifest_f32 = nipizd_manifest(recipe_f32; grid=BoxModelGrid(Float32))
    manifest_f64 = nipizd_manifest(recipe_f32; grid=BoxModelGrid(Float64))
    @test manifest_f32.scalar_type === Float32
    @test manifest_f64.scalar_type === Float64
    @test manifest_f32.sinking_tracers.D isa Float32
    @test manifest_f64.sinking_tracers.D isa Float64

    mktemp() do path, io
        close(io)
        @test export_recipe(path, recipe) == path
        @test import_recipe(path) == recipe
    end

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
    @test occursin("Unsupported Agate recipe schema", sprint(showerror, v1_error))
    @test occursin("agate.model_recipe.v2", sprint(showerror, v1_error))

    invalid_documents = (
        modified(encoded, x -> (x["extra"] = true)),
        modified(encoded, x -> delete!(x["recipe"], "open_bottom")),
        modified(encoded, x -> (x["model"]["family"] = "UnknownModel")),
        modified(encoded, x -> (x["recipe"]["scalar_type"] = "Float32")),
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
