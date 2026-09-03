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

@testset "Recipe provenance does not retain URL credentials" begin
    url = "https://user:token@example.com/org/repo.git?access_token=secret#fragment"
    sanitized = Agate.Construction._sanitize_repository_url(url)
    @test !any(secret -> occursin(secret, sanitized), ("user", "token", "secret"))
end

@testset "NiPiZD versioned family recipes" begin
    family = NiPiZD.NiPiZDFamily()
    @test definition_version(family) == v"0.1.0"
    direct = NiPiZD.construct(; grid=BoxModelGrid(Float32))
    with_recipe, default_recipe = NiPiZD.construct_plus_recipe(; grid=BoxModelGrid(Float32))
    family_constructed = Agate.Construction.construct(family;
        plankton_pfts=default_recipe.plankton_pfts,
        grid=BoxModelGrid(Float32))
    replayed_default = NiPiZD.construct_from_recipe(default_recipe; grid=BoxModelGrid(Float32))
    @test required_biogeochemical_tracers(direct) == required_biogeochemical_tracers(with_recipe)
    @test required_biogeochemical_tracers(replayed_default) == required_biogeochemical_tracers(direct)
    @test required_biogeochemical_tracers(family_constructed) == required_biogeochemical_tracers(direct)
    @test family_constructed.parameters == with_recipe.parameters == direct.parameters == replayed_default.parameters

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
    @test encoded["schema"] == Agate.Construction.recipe_schema() == "agate.model_recipe.v1"
    @test encoded["family"] == "NiPiZD"
    @test encoded["definition_version"] == "0.1.0"
    @test Set(keys(encoded["realization"])) == Set((
        "plankton_pfts",
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
    @test keys(recipe.plankton_pfts) == (:P, :Z)
    @test keys(recipe.plankton_pfts.P) == (:diat,)
    @test keys(recipe.plankton_pfts.Z) == (:microzoo,)
    @test recipe.parameter_overrides == merge(
        inputs.parameters, (palatability_matrix=inputs.palatability_matrix,)
    )
    @test !recipe.open_bottom
    @test recipe.sinking_tracers == inputs.sinking_tracers
    @test decoded == recipe

    mapping_a = (P=(small=[2.0, 1.0], large=[3.0]), Z=(Z=[10.0],))
    mapping_b = (Z=(Z=[10.0],), P=(large=[3.0], small=[1.0, 2.0]))
    overrides_a = (alpha=(small_2=0.3,), maximum_growth_rate=(large_1=1.0e-5, small_1=2.0e-5))
    overrides_b = (maximum_growth_rate=(small_1=2.0e-5, large_1=1.0e-5), alpha=(small_2=0.3,))
    sinking_a = (large_1=0.1, D=2.0)
    sinking_b = (D=2.0, large_1=0.1)
    recipe_a = Agate.Construction.capture_model_recipe(
        family; plankton_pfts=mapping_a, parameter_overrides=overrides_a, sinking_tracers=sinking_a
    )
    recipe_b = Agate.Construction.capture_model_recipe(
        family; plankton_pfts=mapping_b, parameter_overrides=overrides_b, sinking_tracers=sinking_b
    )
    @test recipe_a == recipe_b && isequal(recipe_a, recipe_b)
    @test hash(recipe_a) == hash(recipe_b)
    @test Dict(recipe_a => :same)[recipe_b] === :same
    @test encode_recipe(recipe_a)["content_hash"] == encode_recipe(recipe_b)["content_hash"]

    manifest_a, manifest_b = nipizd_manifest(recipe_a), nipizd_manifest(recipe_b)
    @test manifest_a == manifest_b && isequal(manifest_a, manifest_b)
    @test hash(manifest_a) == hash(manifest_b)
    @test manifest_b in Set([manifest_a])

    direct_multi = Agate.Construction.construct(Agate.ModelDefinition(family); plankton_pfts=mapping_a)
    family_multi = Agate.Construction.construct(family; plankton_pfts=mapping_a)
    @test Agate.Introspection.pfts(direct_multi) == Agate.Introspection.pfts(family_multi)
    @test required_biogeochemical_tracers(direct_multi) == required_biogeochemical_tracers(family_multi)

    @test manifest.pft_entities == (
        diat=(:diat_1, :diat_2),
        microzoo=(:microzoo_1, :microzoo_2),
    )
    @test decoded_manifest == manifest
    @test decoded_manifest.sinking_tracers.D isa Float32

    unsized_recipe = Agate.Construction.ModelRecipe(
        recipe.family,
        recipe.definition_version,
        merge(recipe.plankton_pfts, (P=(diat=nothing,),)),
        recipe.parameter_overrides,
        recipe.sinking_tracers,
        recipe.open_bottom,
    )
    encoded_unsized = encode_recipe(unsized_recipe)
    @test encoded_unsized["realization"]["plankton_pfts"]["P"]["diat"] === nothing
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

    # Identical scientific inputs produce deterministic identity independent of runtime precision.
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
        pop!(x["realization"]["plankton_pfts"])
    end
    invalid_law = rehashed(encoded) do x
        x["realization"]["parameter_overrides"]["linear_mortality"]["law"] = "unknown"
    end
    invalid_matrix = rehashed(encoded) do x
        x["realization"]["parameter_overrides"]["palatability_matrix"][2] = [0.3]
    end
    for document in (invalid_schema, invalid_realization, invalid_law, invalid_matrix)
        @test_throws ArgumentError decode_recipe(document)
    end
    @test_throws ArgumentError decode_recipe("not a recipe document")
end
