using Agate.Construction: decode_recipe, encode_recipe, export_recipe, import_recipe
using Agate.Library.Allometry: AllometricParam, ConstantParam, PowerLaw
using Agate.Models: NiPiZD
using Oceananigans.Units: day
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

    size_structure = (;
        phytoplankton=(diat=Float32[2, 8],),
        zooplankton=(;
            microzoo=(n=2, min_esd=30.0f0, max_esd=90.0f0, splitting=:log_splitting),
        ),
    )
    parameters = (;
        maximum_growth_rate=(diat_2=1.25f0 / day,),
        linear_mortality=AllometricParam(
            PowerLaw(); prefactor=0.05f0 / day, exponent=-0.1f0
        ),
        alpha=ConstantParam(0.2f0 / day),
    )

    _, recipe, manifest = NiPiZD.construct_with_manifest(;
        size_structure,
        scalar_type=Float32,
        parameters,
        palatability_matrix=Float32[0.1 0.9; 0.3 0.7],
        sinking_tracers=(D=2.5f0 / day,),
        open_bottom=false,
    )

    encoded = encode_recipe(recipe)
    decoded = decode_recipe(encoded)
    _, decoded_manifest = NiPiZD.construct_with_manifest(decoded)

    @test encoded["schema"] == "agate.model_recipe.v1"
    @test Set(keys(encoded)) == Set(("schema", "model", "recipe"))
    @test decoded == recipe
    @test decoded_manifest == manifest

    mktemp() do path, io
        close(io)
        @test export_recipe(path, recipe) == path
        @test import_recipe(path) == recipe
    end

    _, nonfinite_recipe = NiPiZD.construct_with_recipe(;
        scalar_type=Float32, palatability_matrix=Float32[NaN 0.9; 0.3 0.7]
    )
    @test_throws ArgumentError encode_recipe(nonfinite_recipe)

    invalid_documents = (
        modified(encoded, x -> (x["extra"] = true)),
        modified(encoded, x -> delete!(x["recipe"], "open_bottom")),
        modified(encoded, x -> (x["schema"] = "agate.model_recipe.v2")),
        modified(encoded, x -> (x["model"]["family"] = "UnknownModel")),
        modified(encoded, x -> (x["recipe"]["scalar_type"] = "Float16")),
        modified(encoded, x -> (
            encoded_named_value(x["recipe"]["parameter_overrides"], :linear_mortality)["law"] = "unknown"
        )),
        modified(encoded, x -> (x["recipe"]["interaction_definitions"][1]["deriver"] = "unknown")),
        modified(encoded, x -> (
            encoded_named_value(x["recipe"]["interaction_overrides"], :palatability_matrix)[2] = [0.3]
        )),
    )

    for document in invalid_documents
        @test_throws ArgumentError decode_recipe(document)
    end
    @test_throws ArgumentError decode_recipe("not a recipe document")
end
