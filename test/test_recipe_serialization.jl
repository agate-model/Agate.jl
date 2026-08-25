using Agate.Construction: decode_recipe, encode_recipe, export_recipe, import_recipe
using Agate.Models: NiPiZD
using Agate.Configuration: Population, Pool
using Agate.Processes:
    AbstractFormulation, Growth, Light, NutrientResponse, Smith, Monod, FrankTNorm,
    formulation_recipe_fields
using OceanBioME: BoxModelGrid
using Test


struct MissingRecipeFields <: AbstractFormulation
    exponent::Float64
end

struct ExplicitRecipeFields <: AbstractFormulation
    exponent::Float64
end

Agate.Processes.formulation_recipe_fields(formulation::ExplicitRecipeFields) =
    (exponent=formulation.exponent,)

function encoded_named_value(encoded, name)
    entry = only(item for item in encoded["entries"] if item["name"] == String(name))
    return entry["value"]
end

function modified(document, mutation::Function)
    copy = deepcopy(document)
    mutation(copy)
    return copy
end

modified(mutation::Function, document) = modified(document, mutation)

function rehash!(document)
    family = Symbol(document["model"]["family"])
    document["recipe_hash"] = Agate.Construction._recipe_hash(family, document["recipe"])
    return document
end

function rehashed(document, mutation)
    copy = modified(document, mutation)
    return rehash!(copy)
end

explicit_json_value(::Nothing) = true
explicit_json_value(::Union{Bool,Int,AbstractFloat,String}) = true
explicit_json_value(x::Dict{String,Any}) = all(explicit_json_value, values(x))
explicit_json_value(x::Vector{Any}) = all(explicit_json_value, x)
explicit_json_value(::Any) = false

@testset "NiPiZD recipe serialization" begin
    _, default_recipe = NiPiZD.construct_plus_recipe()
    default_encoded = encode_recipe(default_recipe)
    @test decode_recipe(default_encoded) == default_recipe
    @test default_encoded["schema"] == "agate.model_recipe.v0.7"

    bgc, recipe = NiPiZD.construct_plus_recipe(; authored_nipizd_inputs(Float32)...)
    manifest = nipizd_manifest(recipe; scalar_type=Float32)

    encoded = encode_recipe(recipe)
    realization = encoded["recipe"]["realization"]
    recipe_hash = encoded["recipe_hash"]
    @test explicit_json_value(encoded)
    bgc.parameters.palatability_matrix[1, 1] = 0f0
    @test recipe.parameter_overrides.palatability_matrix[1, 1] == 0.8f0
    @test encode_recipe(recipe)["recipe_hash"] == recipe_hash
    decoded = decode_recipe(encoded)
    decoded_manifest = nipizd_manifest(decoded; scalar_type=Float32)

    @test Set(keys(encoded)) == Set(("schema", "model", "provenance", "recipe", "recipe_hash"))
    @test Set(keys(encoded["recipe"])) == Set(("components", "processes", "realization"))
    @test Set(keys(realization)) == Set((
        "community",
        "population_groups",
        "parameter_overrides",
        "sinking_tracers",
        "open_bottom",
    ))
    growth_data = encoded["recipe"]["processes"]["growth_P"]
    @test growth_data["factors"] == Dict(
        "light" => Dict(
            "kind" => "light",
            "formulation" => "smith",
            "drivers" => Dict("driver" => "PAR"),
            "bindings" => Dict(
                "alpha" => "alpha",
                "maximum_rate" => "maximum_growth_rate",
            ),
        ),
        "nutrients" => Dict(
            "kind" => "nutrient_response",
            "formulation" => "monod",
            "participants" => Dict("resource" => "N"),
            "bindings" => Dict(
                "K" => Dict("N" => "nutrient_half_saturation"),
            ),
        ),
    )
    phytoplankton_data = encoded["recipe"]["components"]["P"]
    @test phytoplankton_data["kind"] == "population"
    @test phytoplankton_data["states"] == Dict("nitrogen" => "nitrogen")
    @test !haskey(phytoplankton_data, "currency")

    grazing_data = encoded["recipe"]["processes"]["grazing_Z_on_P"]
    @test grazing_data["kind"] == "consumption"
    @test grazing_data["participants"] == Dict("consumer" => "Z", "resource" => "P")
    @test grazing_data["unassimilated_products"] == "D"
    @test grazing_data["bindings"] == Dict(
        "assimilation" => "assimilation_matrix",
        "half_saturation" => "holling_half_saturation",
        "maximum_rate" => "maximum_predation_rate",
        "palatability" => "palatability_matrix",
    )
    mortality_data = encoded["recipe"]["processes"]["linear_mortality_P"]
    @test mortality_data["products"] == Dict(
        "targets" => Dict("detritus" => "D", "nutrient" => "N"),
        "fractions" => Dict("nutrient" => "mortality_export_fraction"),
    )
    @test encoded["recipe"]["processes"]["remineralization_D"]["bindings"] ==
        Dict("rate" => Dict("D" => "detritus_remineralization"))
    @test encoded["provenance"]["agate"]["package"] == "Agate"
    @test encoded["provenance"]["agate"]["version"] == string(Base.pkgversion(Agate))
    @test startswith(recipe_hash, "sha256:")
    @test length(recipe_hash) == 71
    @test decoded == recipe
    @test decoded_manifest == manifest
    @test decoded_manifest.sinking_tracers.D isa Float32

    reordered_growth = Growth(;
        populations=:P,
        factors=(
            nutrients=NutrientResponse(
                Monod(); resource=:N, bindings=(K=:nutrient_half_saturation,)
            ),
            light=Light(
                Smith(); driver=:PAR, bindings=(maximum_rate=:maximum_growth_rate,)
            ),
        ),
    )
    reordered_processes = merge(recipe.processes, (growth_P=reordered_growth,))
    reordered_recipe = Agate.Construction.ProcessModelRecipe(
        recipe.family,
        recipe.components,
        reordered_processes,
        recipe.population_groups,
        recipe.community,
        recipe.parameter_overrides,
        recipe.sinking_tracers,
        recipe.open_bottom,
    )
    @test encode_recipe(reordered_recipe)["recipe_hash"] == recipe_hash

    structured_pool_recipe = Agate.Construction.ProcessModelRecipe(
        recipe.family,
        merge(recipe.components, (POM=Pool(:nitrogen; size_structure=[0.5, 5.0, 50.0]),)),
        recipe.processes,
        recipe.population_groups,
        recipe.community,
        recipe.parameter_overrides,
        recipe.sinking_tracers,
        recipe.open_bottom,
    )
    structured_pool = encode_recipe(structured_pool_recipe)["recipe"]["components"]["POM"]
    @test structured_pool["kind"] == "pool"
    @test structured_pool["size_structure"]["diameters"] == [0.5, 5.0, 50.0]
    @test !haskey(structured_pool, "sinking")

    multistate_recipe = Agate.Construction.ProcessModelRecipe(
        recipe.family,
        merge(
            recipe.components,
            (P=Population(; states=(carbon=:carbon, nitrogen=:nitrogen)),),
        ),
        recipe.processes,
        recipe.population_groups,
        recipe.community,
        recipe.parameter_overrides,
        recipe.sinking_tracers,
        recipe.open_bottom,
    )
    multistate_population = Agate.Construction._component_recipe_data(
        :P, multistate_recipe.components.P, multistate_recipe
    )
    @test multistate_population["states"] ==
        Dict("carbon" => "carbon", "nitrogen" => "nitrogen")

    @test formulation_recipe_fields(FrankTNorm()) == NamedTuple()
    @test_throws ArgumentError Agate.Construction._validated_formulation_recipe_fields(
        MissingRecipeFields(2.0)
    )
    @test Agate.Construction._validated_formulation_recipe_fields(
        ExplicitRecipeFields(2.0)
    ) == (exponent=2.0,)

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

    component_tamper = modified(encoded) do x
        x["recipe"]["components"]["P"]["states"]["nitrogen"] = "carbon"
    end
    @test_throws ArgumentError decode_recipe(component_tamper)

    binding_tamper = modified(encoded) do x
        x["recipe"]["processes"]["growth_P"]["factors"]["light"]["bindings"][
            "maximum_rate"
        ] = "tampered_parameter"
    end
    @test_throws ArgumentError decode_recipe(binding_tamper)

    self_consistent_wrong_contract = rehash!(deepcopy(binding_tamper))
    @test_throws ArgumentError decode_recipe(self_consistent_wrong_contract)

    invalid_schema = modified(encoded, x -> (x["schema"] = "agate.model_recipe.invalid"))
    @test_throws ArgumentError decode_recipe(invalid_schema)

    invalid_documents = (
        modified(encoded, x -> (x["extra"] = true)),
        rehashed(encoded, x -> delete!(x["recipe"]["realization"], "open_bottom")),
        rehashed(encoded, x -> (x["model"]["family"] = "UnknownModel")),
        rehashed(encoded, x -> (x["recipe"]["realization"]["scalar_type"] = "Float32")),
        rehashed(encoded, x -> (
            encoded_named_value(
                x["recipe"]["realization"]["parameter_overrides"], :linear_mortality
            )["law"] = "unknown"
        )),
        rehashed(encoded, x -> (
            encoded_named_value(
                x["recipe"]["realization"]["parameter_overrides"], :palatability_matrix
            )[2] = [0.3]
        )),
    )

    for document in invalid_documents
        @test_throws ArgumentError decode_recipe(document)
    end
    @test_throws ArgumentError decode_recipe("not a recipe document")
end
