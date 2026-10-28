using Agate.Construction: construct, construct_plus_recipe, construct_plus_manifest,
    capture_model_recipe, encode_recipe, decode_recipe, export_recipe, import_recipe
using Agate.Library.Allometry: AllometricParam, SplitPowerLaw
using OceanBioME: BoxModelGrid
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers
using Test

# A model-family fixture defined entirely through Agate's public extension interface.
struct TestNPZFamily <: Agate.ModelFamilies.AbstractModelFamily end
Agate.Construction.family_id(::TestNPZFamily) = :TestNPZ
Agate.ModelFamilies.default_components(::TestNPZFamily) = runtime_npz_definition().components
Agate.ModelFamilies.default_processes(::TestNPZFamily) = runtime_npz_definition().processes
Agate.Parameters.parameter_definitions(::TestNPZFamily) = runtime_npz_definition().parameters

struct AlternativeTestNPZFamily <: Agate.ModelFamilies.AbstractModelFamily end
Agate.Construction.family_id(::AlternativeTestNPZFamily) = :TestNPZ
Agate.ModelFamilies.default_components(::AlternativeTestNPZFamily) = runtime_npz_definition().components
Agate.ModelFamilies.default_processes(::AlternativeTestNPZFamily) = runtime_npz_definition().processes
Agate.Parameters.parameter_definitions(::AlternativeTestNPZFamily) = merge(
    runtime_npz_definition().parameters,
    (detritus_remineralization=Agate.Parameters.Parameter(2e-5),),
)

modified(f::Function, x) = (copy=deepcopy(x); f(copy); copy)
function rehashed(f::Function, x)
    return modified(x) do document
        f(document)
        document["content_hash"] = Agate.Construction._recipe_hash(
            Symbol(document["family"]), document["realization"]
        )
    end
end

@testset "Generic model-family recipe serialization and replay" begin
    family = TestNPZFamily()
    pfts = (P=(P=[2.0, 10.0],), Z=(Z=[20.0, 100.0],))
    grid = BoxModelGrid(Float32)
    bgc, default_recipe = construct_plus_recipe(family; plankton_pfts=pfts, grid)
    replayed = construct(family, default_recipe; grid)
    direct = construct(runtime_npz_definition(); plankton_pfts=pfts, grid)
    @test bgc.parameters == replayed.parameters == direct.parameters
    @test required_biogeochemical_tracers(bgc) == required_biogeochemical_tracers(replayed)

    # Provider selection supplies scientific defaults; the recipe supplies configuration.
    other = construct(AlternativeTestNPZFamily(), default_recipe; grid)
    @test other.parameters.detritus_remineralization == Float32(2e-5)
    @test other.parameters.detritus_remineralization != bgc.parameters.detritus_remineralization

    pfts = (
        P=(diat=[2.0, 8.0],),
        Z=(microzoo=[30.0, 90.0],),
    )
    inputs = (;
        plankton_pfts=pfts,
        parameter_overrides=(;
            maximum_growth_rate=(diat_2=1.25e-5,),
            palatability_matrix=Float32[0.8 0.2; 0.3 0.7],
        ),
        sinking_tracers=(D=2.5 / 86400,),
        open_bottom=false,
    )
    original, recipe = construct_plus_recipe(family; inputs..., grid)
    encoded = encode_recipe(recipe)
    decoded = decode_recipe(encoded)
    @test recipe.family == :TestNPZ
    @test decoded == recipe
    @test encoded["schema"] == Agate.Construction.recipe_schema() == "agate.model_recipe.v0.3"
    @test encoded["family"] == "TestNPZ"
    @test Set(keys(encoded)) == Set(("schema", "family", "realization", "content_hash"))
    @test Set(keys(encoded["realization"])) == Set((
        "plankton_pfts", "parameter_overrides", "setting_overrides", "sinking_tracers", "open_bottom",
    ))
    @test startswith(encoded["content_hash"], "sha256:")
    @test length(encoded["content_hash"]) == 71
    @test keys(recipe.plankton_pfts) == (:P, :Z)
    @test recipe.parameter_overrides == inputs.parameter_overrides
    @test recipe.sinking_tracers == inputs.sinking_tracers
    @test !recipe.open_bottom

    reconstructed, manifest = construct_plus_manifest(family, decoded; grid)
    @test reconstructed.parameters == original.parameters
    @test manifest.tracer_order == required_biogeochemical_tracers(original)
    @test manifest.pft_entities == (diat=(:diat_1, :diat_2), microzoo=(:microzoo_1, :microzoo_2))
    @test manifest.scalar_type == Float32
    @test manifest.sinking_tracers.D isa Float32

    mktemp() do path, io
        close(io)
        @test export_recipe(path, recipe) == path
        @test import_recipe(path) == recipe
        @test construct(family, import_recipe(path); grid).parameters == original.parameters
    end

    # Canonical order and content fingerprint describe configuration, not runtime precision.
    pfts_a = (P=(small=[2.0, 1.0], large=[3.0]), Z=(Z=[10.0],))
    pfts_b = (Z=(Z=[10.0],), P=(large=[3.0], small=[1.0, 2.0]))
    override_a = (maximum_growth_rate=(small_2=2e-5, large_1=1e-5),)
    override_b = (maximum_growth_rate=(large_1=1e-5, small_2=2e-5),)
    recipe_a = capture_model_recipe(family;
        plankton_pfts=pfts_a, parameter_overrides=override_a, sinking_tracers=(D=2.0,))
    recipe_b = capture_model_recipe(family;
        plankton_pfts=pfts_b, parameter_overrides=override_b, sinking_tracers=(D=2.0,))
    @test recipe_a == recipe_b && isequal(recipe_a, recipe_b)
    @test Dict(recipe_a => :same)[recipe_b] == :same
    @test encode_recipe(recipe_a)["content_hash"] == encode_recipe(recipe_b)["content_hash"]
    @test encode_recipe(decode_recipe(encode_recipe(recipe_a))) == encode_recipe(recipe_a)

    allometric = capture_model_recipe(family;
        plankton_pfts=pfts,
        parameter_overrides=(maximum_growth_rate=AllometricParam(
            SplitPowerLaw(); prefactor=1.2e-5, breakpoint=3.0,
            small_exponent=0.28, large_exponent=-0.15,
        ),),
    )
    @test decode_recipe(encode_recipe(allometric)) == allometric

    _, float32_recipe = construct_plus_recipe(family; plankton_pfts=pfts, grid=BoxModelGrid(Float32))
    _, float64_recipe = construct_plus_recipe(family; plankton_pfts=pfts, grid=BoxModelGrid(Float64))
    @test encode_recipe(float32_recipe) == encode_recipe(float64_recipe)
    _, manifest64 = construct_plus_manifest(family, float32_recipe; grid=BoxModelGrid(Float64))
    @test manifest64.scalar_type == Float64

    # Decode is independent of any global family registry; the chosen provider checks identity.
    wrong_family = rehashed(encoded) do document
        document["family"] = "OtherModel"
    end
    @test decode_recipe(wrong_family).family == :OtherModel
    @test_throws r"Recipe family" construct(family, decode_recipe(wrong_family))
    @test_throws ArgumentError decode_recipe("not a recipe document")
    @test_throws ArgumentError decode_recipe(modified(encoded) do document
        document["realization"]["open_bottom"] = true
    end)

    invalid_schema = modified(encoded) do document
        document["schema"] = "agate.model_recipe.v0.2"
    end
    @test_throws r"Unsupported Agate recipe schema" decode_recipe(invalid_schema)
    invalid_diameter = rehashed(encoded) do document
        document["realization"]["plankton_pfts"]["P"]["diat"]["diameters"][1] = 0.0
    end
    @test_throws ArgumentError decode_recipe(invalid_diameter)

    # Invalid content fails atomically, leaving the existing recipe file unchanged.
    mktemp() do path, io
        write(io, "original")
        close(io)
        invalid = Agate.Construction.ModelRecipe(
            recipe.family, recipe.plankton_pfts, (unsupported=Dict(:key => 1),),
            recipe.setting_overrides, recipe.sinking_tracers, recipe.open_bottom,
        )
        @test_throws ArgumentError export_recipe(path, invalid)
        @test read(path, String) == "original"
    end
end
