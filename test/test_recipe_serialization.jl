using Agate.Construction: decode_recipe, encode_recipe, export_recipe, import_recipe
using Agate.ModelFamilies: definition_version
using Agate.Models: NiPiZD
using OceanBioME: BoxModelGrid
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers, biogeochemical_drift_velocity
using Test

function modified(mutation::Function, document)
    copy = deepcopy(document)
    mutation(copy)
    return copy
end

function rehash!(document)
    family = Symbol(document["family"])
    version = VersionNumber(document["definition_version"])
    document["content_hash"] = Agate.Construction._recipe_hash(
        family, version, document["realization"]
    )
    return document
end

rehashed(mutation::Function, document) = rehash!(modified(mutation, document))

explicit_json_value(::Nothing) = true
explicit_json_value(::Union{Bool,Int,AbstractFloat,String}) = true
explicit_json_value(x::Dict{String,Any}) = all(explicit_json_value, values(x))
explicit_json_value(x::Vector{Any}) = all(explicit_json_value, x)
explicit_json_value(::Any) = false

@testset "NiPiZD versioned family recipes" begin
    family = NiPiZD.NiPiZDFamily()
    @test definition_version(family) == v"0.1.0"

    direct = NiPiZD.construct(; grid=BoxModelGrid(Float32))
    with_recipe, default_recipe = NiPiZD.construct_plus_recipe(; grid=BoxModelGrid(Float32))
    replayed_default = NiPiZD.construct_from_recipe(default_recipe; grid=BoxModelGrid(Float32))
    @test required_biogeochemical_tracers(direct) == required_biogeochemical_tracers(with_recipe)
    @test required_biogeochemical_tracers(replayed_default) == required_biogeochemical_tracers(direct)
    @test with_recipe.parameters == direct.parameters == replayed_default.parameters

    inputs = authored_nipizd_inputs(Float32)
    bgc, recipe = NiPiZD.construct_plus_recipe(; inputs...)
    manifest = nipizd_manifest(recipe; scalar_type=Float32)
    encoded = encode_recipe(recipe)
    decoded = decode_recipe(encoded)
    decoded_manifest = nipizd_manifest(decoded; scalar_type=Float32)

    @test Set(keys(encoded)) == Set((
        "schema",
        "family",
        "definition_version",
        "realization",
        "provenance",
        "content_hash",
    ))
    @test encoded["schema"] == "agate.model_recipe.v1"
    @test encoded["family"] == "NiPiZD"
    @test encoded["definition_version"] == "0.1.0"
    @test Set(keys(encoded["realization"])) == Set((
        "plankton_pfts",
        "pft_size_structures",
        "parameter_overrides",
        "sinking_tracers",
        "open_bottom",
    ))
    @test !haskey(encoded, "model")
    @test !haskey(encoded, "recipe")
    @test explicit_json_value(encoded)
    @test startswith(encoded["content_hash"], "sha256:")
    @test length(encoded["content_hash"]) == 71
    @test encoded["provenance"]["agate"]["package"] == "Agate"
    @test encoded["provenance"]["agate"]["version"] == string(Base.pkgversion(Agate))

    @test recipe.family === :NiPiZD
    @test recipe.definition_version == definition_version(family)
    @test recipe.plankton_pfts == (P=(:diat,), Z=(:microzoo,))
    @test keys(recipe.pft_size_structures) == (:diat, :microzoo)
    @test recipe.parameter_overrides == merge(
        inputs.parameters, (palatability_matrix=inputs.palatability_matrix,)
    )
    @test !recipe.open_bottom
    @test recipe.sinking_tracers == inputs.sinking_tracers
    @test decoded == recipe
    @test manifest.pft_entities == (
        diat=(:diat_1, :diat_2),
        microzoo=(:microzoo_1, :microzoo_2),
    )
    @test decoded_manifest == manifest
    @test decoded_manifest.sinking_tracers.D isa Float32

    unsized_recipe = Agate.Construction.ModelRecipe(
        recipe.family,
        recipe.definition_version,
        recipe.plankton_pfts,
        (diat=nothing, microzoo=recipe.pft_size_structures.microzoo),
        recipe.parameter_overrides,
        recipe.sinking_tracers,
        recipe.open_bottom,
    )
    encoded_unsized = encode_recipe(unsized_recipe)
    diat_size = only(
        entry for entry in encoded_unsized["realization"]["pft_size_structures"]
        if entry["pft"] == "diat"
    )
    @test diat_size["diameters"] === nothing
    @test decode_recipe(encoded_unsized) == unsized_recipe

    # Captured inputs are isolated from the constructed runtime model.
    recipe_hash = encoded["content_hash"]
    bgc.parameters.palatability_matrix[1, 1] = 0f0
    @test recipe.parameter_overrides.palatability_matrix[1, 1] == 0.8f0
    @test encode_recipe(recipe)["content_hash"] == recipe_hash

    replayed = NiPiZD.construct_from_recipe(decoded; grid=BoxModelGrid(Float32))
    @test all(
        getproperty(replayed.parameters, name) == getproperty(decoded_manifest.parameters, name)
        for name in keys(replayed.parameters)
    )
    @test !hasproperty(replayed.parameters, :specificity)
    @test hasproperty(decoded_manifest.parameters, :specificity)
    @test required_biogeochemical_tracers(replayed) == decoded_manifest.tracer_order
    original_drift = biogeochemical_drift_velocity(bgc, Val(:D))
    replayed_drift = biogeochemical_drift_velocity(replayed, Val(:D))
    @test replayed_drift.w.data == original_drift.w.data

    # Identical ordered inputs produce deterministic scientific identity independent of runtime precision.
    _, recipe_f32 = NiPiZD.construct_plus_recipe(;
        grid=BoxModelGrid(Float32), sinking_tracers=(D=2.5 / 86400,), open_bottom=false
    )
    _, recipe_f64 = NiPiZD.construct_plus_recipe(;
        grid=BoxModelGrid(Float64), sinking_tracers=(D=2.5 / 86400,), open_bottom=false
    )
    encoded_f32 = encode_recipe(recipe_f32)
    encoded_f64 = encode_recipe(recipe_f64)
    @test recipe_f32 == recipe_f64
    @test encoded_f32 == encoded_f64

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

    version_warning = modified(encoded) do x
        x["provenance"]["agate"]["version"] = "0.9.0"
    end
    @test_logs (:warn, r"Agate version differs") decode_recipe(version_warning)

    stale_hash = modified(encoded) do x
        x["realization"]["open_bottom"] = !x["realization"]["open_bottom"]
    end
    @test_throws ArgumentError decode_recipe(stale_hash)

    unknown_family = rehashed(encoded) do x
        x["family"] = "UnknownModel"
    end
    @test_throws ArgumentError decode_recipe(unknown_family)

    missing_version = modified(encoded) do x
        delete!(x, "definition_version")
    end
    @test_throws ArgumentError decode_recipe(missing_version)

    mismatched_version = rehashed(encoded) do x
        x["definition_version"] = "0.1.1"
    end
    @test_throws ArgumentError decode_recipe(mismatched_version)

    bumped_recipe = Agate.Construction.ModelRecipe(
        recipe.family,
        v"0.1.1",
        recipe.plankton_pfts,
        recipe.pft_size_structures,
        recipe.parameter_overrides,
        recipe.sinking_tracers,
        recipe.open_bottom,
    )
    @test encode_recipe(bumped_recipe)["content_hash"] != encoded["content_hash"]
    @test_throws ArgumentError NiPiZD.construct_from_recipe(bumped_recipe)
    @test required_biogeochemical_tracers(NiPiZD.construct()) ==
          (:N, :D, :P_1, :P_2, :Z_1, :Z_2)

    invalid_schema = modified(encoded) do x
        x["schema"] = "agate.model_recipe.invalid"
    end
    invalid_realization = rehashed(encoded) do x
        pop!(x["realization"]["pft_size_structures"])
    end
    invalid_law = rehashed(encoded) do x
        override = only(
            entry for entry in x["realization"]["parameter_overrides"]
            if entry["name"] == "linear_mortality"
        )
        override["value"]["law"] = "unknown"
    end
    invalid_matrix = rehashed(encoded) do x
        override = only(
            entry for entry in x["realization"]["parameter_overrides"]
            if entry["name"] == "palatability_matrix"
        )
        override["value"][2] = [0.3]
    end
    for document in (invalid_schema, invalid_realization, invalid_law, invalid_matrix)
        @test_throws ArgumentError decode_recipe(document)
    end
    @test_throws ArgumentError decode_recipe("not a recipe document")
end
