using Agate.Construction: decode_recipe, encode_recipe, export_recipe, import_recipe
using Agate.Library.Allometry: AllometricParam, ConstantParam, PowerLaw
using Agate.Models: NiPiZD
using Oceananigans.Units: day
using Test

struct InvalidIdentifierDeriver <: Agate.Configuration.AbstractMatrixDeriver end
Agate.Configuration.matrix_deriver_identifier(::InvalidIdentifierDeriver) = "invalid"

struct StructuralEqualityFixture
    values::Vector{Int}
end

function recipe_with_interaction_definitions(recipe, interaction_definitions)
    return Agate.Construction.ModelRecipe(
        recipe.family,
        recipe.community,
        recipe.parameter_definitions,
        recipe.parameter_overrides,
        interaction_definitions,
        recipe.interaction_overrides,
        recipe.ecological_roles,
        recipe.interaction_roles,
        recipe.parameter_roles,
        recipe.auxiliary_fields,
        recipe.sinking_tracers,
        recipe.open_bottom,
        recipe.scalar_type,
    )
end

function encoded_named_entry(encoded, name)
    return only(filter(item -> item["name"] == String(name), encoded["entries"]))
end

encoded_named_value(encoded, name) = encoded_named_entry(encoded, name)["value"]

@testset "NiPiZD recipe serialization" begin
    _, default_recipe = NiPiZD.construct_with_recipe()
    default_encoded = encode_recipe(default_recipe)
    @test default_recipe.family === :NiPiZD
    @test default_encoded["recipe"]["interaction_overrides"] ==
          Dict{String,Any}("type" => "named_tuple", "entries" => Any[])
    @test decode_recipe(default_encoded) == default_recipe

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
    @test encoded["schema"] == "agate.model_recipe.v1"
    encoded_matrix = encoded_named_value(
        encoded["recipe"]["interaction_overrides"], :palatability_matrix
    )
    @test encoded_matrix isa Vector
    @test all(row -> row isa Vector, encoded_matrix)
    @test all(value -> value isa Float64, Iterators.flatten(encoded_matrix))
    @test Set(keys(encoded)) == Set(("schema", "model", "recipe"))
    @test all(
        !haskey(definition["spec"], "doc")
        for definition in encoded["recipe"]["parameter_definitions"]
    )
    decoded = decode_recipe(encoded)
    _, decoded_manifest = NiPiZD.construct_with_manifest(decoded)
    @test decoded == recipe
    @test decoded_manifest == manifest

    _, nonfinite_recipe = NiPiZD.construct_with_recipe(;
        scalar_type=Float32, palatability_matrix=Float32[NaN 0.9; 0.3 0.7]
    )
    @test_throws ArgumentError encode_recipe(nonfinite_recipe)

    definition = only(
        filter(d -> d.spec.name == :maximum_growth_rate, decoded.parameter_definitions)
    )
    @test definition.spec.materialization.role === :producers
    @test definition.spec.materialization.fill_value == 0

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
    invalid["model"]["family"] = "UnknownModel"
    @test_throws ArgumentError decode_recipe(invalid)

    invalid = deepcopy(encoded)
    invalid["recipe"]["scalar_type"] = "Float16"
    @test_throws ArgumentError decode_recipe(invalid)

    invalid = deepcopy(encoded)
    microzoo = only(filter(item -> item["name"] == "microzoo", invalid["recipe"]["community"]))
    encoded_named_value(microzoo["spec"], :diameters)["splitting"] = "geometric"
    @test_throws ArgumentError decode_recipe(invalid)

    invalid = deepcopy(encoded)
    encoded_named_value(
        invalid["recipe"]["parameter_overrides"], :linear_mortality
    )["relationship"] = "unknown"
    @test_throws ArgumentError decode_recipe(invalid)

    invalid = deepcopy(encoded)
    encoded_named_entry(
        invalid["recipe"]["interaction_overrides"], :palatability_matrix
    )["value"] = 1
    @test_throws ArgumentError decode_recipe(invalid)

    invalid = deepcopy(encoded)
    encoded_named_entry(
        invalid["recipe"]["interaction_overrides"], :palatability_matrix
    )["value"] = [1, 2]
    @test_throws ArgumentError decode_recipe(invalid)

    invalid = deepcopy(encoded)
    encoded_named_entry(
        invalid["recipe"]["interaction_overrides"], :palatability_matrix
    )["value"] = [[0.8, 0.2], [0.3]]
    @test_throws ArgumentError decode_recipe(invalid)

    invalid = deepcopy(encoded)
    invalid["recipe"]["interaction_definitions"][1]["deriver"] = "unknown"
    @test_throws ArgumentError decode_recipe(invalid)

    invalid_identifier_recipe = recipe_with_interaction_definitions(
        default_recipe,
        (invalid=Agate.Configuration.MatrixDefinition(InvalidIdentifierDeriver()),),
    )
    @test_throws ArgumentError encode_recipe(invalid_identifier_recipe)

    manifest_a = Agate.Construction.ModelManifest(
        (fixture=StructuralEqualityFixture([1, 2]),),
        (;),
        (),
        (),
        (),
        (;),
        (consumers=(), prey=()),
        (;),
        (;),
        nothing,
        true,
        Float64,
    )
    manifest_b = Agate.Construction.ModelManifest(
        (fixture=StructuralEqualityFixture([1, 2]),),
        (;),
        (),
        (),
        (),
        (;),
        (consumers=(), prey=()),
        (;),
        (;),
        nothing,
        true,
        Float64,
    )
    @test manifest_a == manifest_b
    manifest_b.parameters.fixture.values[2] = 3
    @test manifest_a != manifest_b
    @test_throws ArgumentError decode_recipe("not a recipe document")
end
