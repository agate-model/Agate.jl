using Agate.Construction: ModelRecipe, capture_model_recipe, decode_recipe, encode_recipe,
                          export_recipe, import_recipe, replay_family
using Agate.ModelFamilies: AbstractModelFamily, definition_version, model_lineage
using Agate.Models: NiPiZD
using Test

struct LineageFixture{L} <: AbstractModelFamily
    lineage::L
end

Agate.Construction.family_id(::LineageFixture) = :LineageFixture
Agate.ModelFamilies.definition_version(::LineageFixture) = v"0.3.0"
Agate.ModelFamilies.model_lineage(family::LineageFixture) = family.lineage
Agate.ModelFamilies.default_components(::LineageFixture) =
    Agate.ModelFamilies.default_components(NiPiZD.NiPiZDFamily())

const TEST_LINEAGE = (:Parent => v"0.2.0", :Original => v"0.1.0")
Agate.Construction.registered_family(::Val{:LineageFixture}) = LineageFixture(TEST_LINEAGE)

function fixture_recipe(lineage=TEST_LINEAGE)
    return capture_model_recipe(
        LineageFixture(lineage);
        plankton_pfts=(P=(diat=nothing,), Z=(microzoo=nothing,)),
    )
end

@testset "Model lineage" begin
    @test isempty(model_lineage(NiPiZD.NiPiZDFamily()))
    recipe = fixture_recipe()
    @test recipe.lineage == TEST_LINEAGE
    document = encode_recipe(recipe)
    @test document["lineage"] == [
        Dict("family" => "Parent", "definition_version" => "0.2.0"),
        Dict("family" => "Original", "definition_version" => "0.1.0"),
    ]
    @test decode_recipe(document) == recipe
    @test replay_family(recipe) isa LineageFixture
    mktemp() do path, io
        close(io)
        export_recipe(path, recipe)
        @test import_recipe(path) == recipe
    end

    changed = deepcopy(document)
    changed["lineage"][1]["family"] = "Other"
    @test_throws "content_hash does not match" decode_recipe(changed)
    changed["content_hash"] = Agate.Construction._recipe_hash(
        recipe.family, recipe.definition_version, changed["lineage"], changed["realization"]
    )
    @test_throws "lineage" decode_recipe(changed)

    other = ModelRecipe(recipe.family, recipe.definition_version, (), recipe.plankton_pfts,
                        recipe.parameter_overrides, recipe.setting_overrides,
                        recipe.sinking_tracers, recipe.open_bottom)
    @test encode_recipe(other)["content_hash"] != document["content_hash"]
    @test_throws "lineage" replay_family(other)

    for lineage in ((:Parent => v"0.2.0", :Parent => v"0.2.0"),
                    (:LineageFixture => v"0.3.0",),
                    ((family=:Parent, definition_version=v"0.2.0"),))
        @test_throws ArgumentError fixture_recipe(lineage)
    end
end
