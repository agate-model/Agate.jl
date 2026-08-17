using Agate.Construction: decode_recipe, encode_recipe, export_recipe, import_recipe
using Agate.Library.Allometry: AllometricParam, ConstantParam, PowerLaw
using Agate.Models: NiPiZD
using Oceananigans.Units: day
using Test

function encoded_named_value(encoded, name)
    entry = only(filter(item -> item["name"] == String(name), encoded["entries"]))
    return entry["value"]
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

    _, recipe, realization = NiPiZD.construct_with_realization(;
        size_structure,
        scalar_type=Float32,
        parameters,
        palatability_matrix=Float32[NaN Inf; -Inf 0.7],
        sinking_tracers=(D=2.5f0 / day,),
        open_bottom=false,
    )

    encoded = encode_recipe(recipe)
    decoded = decode_recipe(encoded)
    _, decoded_realization = NiPiZD.construct_with_realization(decoded)
    @test decoded == recipe
    @test decoded_realization == realization

    mktemp() do path, io
        close(io)
        @test export_recipe(path, recipe) == path
        @test import_recipe(path) == recipe
    end

    invalid = deepcopy(encoded)
    invalid["extra"] = true
    @test_throws ArgumentError decode_recipe(invalid)

    invalid = deepcopy(encoded)
    invalid["schema"] = "agate.model_recipe.v2"
    @test_throws ArgumentError decode_recipe(invalid)

    invalid = deepcopy(encoded)
    invalid["recipe"]["scalar_type"] = "Float16"
    @test_throws ArgumentError decode_recipe(invalid)

    invalid = deepcopy(encoded)
    microzoo = only(filter(item -> item["name"] == "microzoo", invalid["recipe"]["community"]))
    encoded_named_value(microzoo["spec"], :diameters)["splitting"] = "geometric"
    @test_throws ArgumentError decode_recipe(invalid)

    invalid = deepcopy(encoded)
    encoded_named_value(invalid["recipe"]["parameter_overrides"], :linear_mortality)["relationship"] = "unknown"
    @test_throws ArgumentError decode_recipe(invalid)
end
