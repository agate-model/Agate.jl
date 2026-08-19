using Agate
using JSON
using SHA
using Test

using Agate.Construction: encode_recipe
using Agate.Factories: parameter_directory
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers

const CYCLE01_REFERENCE_DIR = joinpath(@__DIR__, "reference")

load_cycle01_reference(name) = JSON.parsefile(joinpath(CYCLE01_REFERENCE_DIR, name))

@testset "v0.11 baseline contract" begin
    baseline = load_cycle01_reference("nipizd_v011_baseline.json")
    trajectory = baseline["numerical"]["historical_npzd_trajectory"]
    trajectory_path = joinpath(CYCLE01_REFERENCE_DIR, trajectory["file"])

    @test baseline["package"] == Dict("name" => "Agate", "version" => "0.11.0")
    @test bytes2hex(sha256(read(trajectory_path))) == trajectory["sha256"]

    rows = filter(
        row -> !isempty(row) && !startswith(row, '#') && row != "day,P",
        readlines(trajectory_path),
    )
    @test length(rows) == trajectory["samples"]

    bgc = Agate.Models.NiPiZD.construct(; grid=dummy_grid(Float64))
    tracers = required_biogeochemical_tracers(bgc)
    @test String.(tracers) == baseline["runtime"]["tracer_order"]
    @test !any(type -> type === Any, fieldtypes(typeof(bgc.tracer_functions)))

    _, recipe = Agate.Models.NiPiZD.construct_plus_recipe()
    @test encode_recipe(recipe)["schema"] == baseline["runtime"]["recipe_schema"]
end

@testset "v0.12 normalization contracts" begin
    contracts = load_cycle01_reference("v012_normalization_contracts.json")
    target = load_cycle01_reference("nipizd_model_recipe_v3_target.json")
    recipe = target["recipe"]
    components = recipe["components"]
    processes = recipe["processes"]
    bindings = recipe["parameter_bindings"]

    @test target["schema"] == contracts["recipe"]["schema"] == "agate.model_recipe.v3"
    @test contracts["processes"]["identity_source"] == "named_process_key"
    @test contracts["processes"]["identity_is_order_independent"]
    @test contracts["processes"]["rate_topology_is_explicit"]
    @test contracts["parameters"]["applicability_source"] == "process_participation"
    @test contracts["stoichiometry"]["ratio_orientation"] == ["from", "to"]
    @test contracts["runtime"]["tendency_execution"] == "compiled_equation_per_tracer"

    @test Set(keys(components)) == Set(("P", "Z", "N", "D"))
    @test all(component -> component["currency"] == "nitrogen", values(components))
    @test components["P"]["kind"] == components["Z"]["kind"] == "population"
    @test components["N"]["kind"] == components["D"]["kind"] == "pool"

    expected_rate_axes = Dict(
        "growth_P" => ["population"],
        "grazing_Z_on_P" => ["consumer", "resource"],
        "linear_mortality_P" => ["population"],
        "linear_mortality_Z" => ["population"],
        "quadratic_mortality_Z" => ["population"],
        "remineralization_D" => ["source"],
    )
    @test Set(keys(processes)) == Set(keys(expected_rate_axes))
    @test all(id -> processes[id]["rate_axes"] == expected_rate_axes[id], keys(processes))

    component_ids = Set(keys(components))
    for process in values(processes), participant in values(process["participants"])
        @test participant in component_ids
    end
    limitation = processes["growth_P"]["subformulations"]["limitation"]
    @test limitation["participants"]["resource"] in component_ids
    @test processes["growth_P"]["drivers"] == Dict("light" => "PAR")

    parameter_names = Set(
        String(spec.name) for spec in parameter_directory(Agate.Models.NiPiZD.NiPiZDFactory())
    )
    @test Set(keys(bindings)) == parameter_names

    requirement_fields = Set(contracts["parameters"]["requirement_identity_fields"])
    for binding in values(bindings), requirement in binding["provides"]
        @test issubset(requirement_fields, Set(keys(requirement)))
        @test haskey(processes, requirement["process"])
        @test requirement["path"] isa Vector
        @test !isempty(requirement["formulation"])
        @test !isempty(requirement["slot"])
    end

    @test bindings["nutrient_half_saturation"]["provides"] == [
        Dict(
            "process" => "growth_P",
            "path" => ["limitation"],
            "formulation" => "monod",
            "slot" => "K",
            "qualifier" => Dict("resource" => "N"),
            "axes" => ["population"],
        ),
    ]
    @test only(bindings["protection"]["provides"])["axes"] == ["resource"]
end
