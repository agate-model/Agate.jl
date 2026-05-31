using Adapt
using JSON
using Oceananigans.Units: day
using OceanBioME: BoxModelGrid
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers

struct UnsupportedManifestValue end


@testset "Manifest value serialization" begin
    value = Agate.Construction.manifest_parameter_value

    @test value(:foo) == "foo"
    @test value((1, :x, NaN)) == Any[1, "x", "NaN"]
    @test value((a=1, b=:x)) == Dict{String,Any}("a" => 1, "b" => "x")
    @test value(Dict(:a => [1, Inf])) == Dict{String,Any}("a" => Any[1, "Inf"])
    @test_throws ArgumentError value(UnsupportedManifestValue())
end

@testset "Model-level construction manifests" begin
    darwin_bgc = Agate.Models.DARWIN.construct(
        ;
        grid=dummy_grid(Float32),
        phyto_size_structure=(n=3, min_esd=1.5, max_esd=20.0, splitting=:log_splitting),
        zoo_size_structure=(n=2, min_esd=20.0, max_esd=100.0, splitting=:log_splitting),
    )
    traced_darwin, darwin_manifest = Agate.Models.DARWIN.construct_with_manifest(
        ;
        grid=dummy_grid(Float32),
        phyto_size_structure=(n=3, min_esd=1.5, max_esd=20.0, splitting=:log_splitting),
        zoo_size_structure=(n=2, min_esd=20.0, max_esd=100.0, splitting=:log_splitting),
    )

    @test typeof(traced_darwin) == typeof(darwin_bgc)
    @test darwin_manifest["model"]["id"] == "DARWIN/default"
    @test darwin_manifest["model"]["family"] == "DARWIN"
    @test !haskey(darwin_manifest, "replayable")
    @test darwin_manifest["resolved"]["scalar_type"] == "Float32"
    @test darwin_manifest["resolved"]["sinking"] == Dict{String,Any}(
        "enabled" => false, "tracers" => nothing, "open_bottom" => true
    )
    @test "P3" in darwin_manifest["resolved"]["tracers"]
    @test darwin_manifest["recipe"]["type"] == "model_constructor"
    @test darwin_manifest["recipe"]["family"] == "DARWIN"
    @test !haskey(darwin_manifest["recipe"], "constructor")
    @test darwin_manifest["recipe"]["kwargs"]["phyto_size_structure"] ==
        darwin_manifest["resolved"]["plankton_diameters_by_group"]["P"]
    @test darwin_manifest["recipe"]["kwargs"]["zoo_size_structure"] ==
        darwin_manifest["resolved"]["plankton_diameters_by_group"]["Z"]
    @test darwin_manifest["recipe"]["kwargs"]["parameters"] isa Dict
    @test darwin_manifest["recipe"]["kwargs"]["sinking_tracers"] === nothing
    @test darwin_manifest["recipe"]["kwargs"]["open_bottom"] == true
    @test darwin_manifest["recipe"]["kwargs"]["scalar_type"] == "Float32"

    _, sinking_manifest = Agate.Models.NiPiZD.construct_with_manifest(
        ;
        grid=BoxModelGrid(),
        sinking_tracers=(P1=0.2551 / day, D=2.7489 / day),
        open_bottom=false,
    )
    @test sinking_manifest["resolved"]["has_sinking_velocities"]
    @test sinking_manifest["resolved"]["sinking"]["enabled"]
    @test sinking_manifest["resolved"]["sinking"]["open_bottom"] == false
    @test sinking_manifest["resolved"]["sinking"]["tracers"] == Dict{String,Any}(
        "P1" => 0.2551 / day, "D" => 2.7489 / day
    )
    @test sinking_manifest["recipe"]["kwargs"]["sinking_tracers"] == Any[
        Dict{String,Any}("name" => "P1", "value" => 0.2551 / day),
        Dict{String,Any}("name" => "D", "value" => 2.7489 / day),
    ]
    @test sinking_manifest["recipe"]["kwargs"]["open_bottom"] == false

    nipizd_bgc = Agate.Models.NiPiZD.construct(; grid=dummy_grid(Float32))
    traced_nipizd, nipizd_manifest = Agate.Models.NiPiZD.construct_with_manifest(
        ; grid=dummy_grid(Float32)
    )

    @test typeof(traced_nipizd) == typeof(nipizd_bgc)
    @test nipizd_manifest["model"]["id"] == "NiPiZD/default"
    @test nipizd_manifest["model"]["family"] == "NiPiZD"
    @test !haskey(nipizd_manifest, "replayable")
    @test nipizd_manifest["resolved"]["scalar_type"] == "Float32"
    @test haskey(nipizd_manifest["resolved"], "parameters")
    @test nipizd_manifest["recipe"]["type"] == "model_constructor"
    @test nipizd_manifest["recipe"]["family"] == "NiPiZD"
    @test !haskey(nipizd_manifest["recipe"], "constructor")
    @test haskey(nipizd_manifest["recipe"]["kwargs"], "parameters")

    replayed_darwin = Agate.Models.construct_from_manifest(darwin_manifest; grid=dummy_grid(Float32))
    @test typeof(replayed_darwin) == typeof(darwin_bgc)
    @test required_biogeochemical_tracers(replayed_darwin) == required_biogeochemical_tracers(darwin_bgc)

    replayed_nipizd = Agate.Models.construct_from_manifest(nipizd_manifest; grid=dummy_grid(Float32))
    @test typeof(replayed_nipizd) == typeof(nipizd_bgc)

    @test required_biogeochemical_tracers(replayed_nipizd) == required_biogeochemical_tracers(nipizd_bgc)

    replayed_sinking_nipizd = Agate.Models.construct_from_manifest(sinking_manifest; grid=BoxModelGrid())
    @test typeof(replayed_sinking_nipizd) == typeof(Agate.Models.NiPiZD.construct(; grid=BoxModelGrid(), sinking_tracers=(P1=0.2551 / day, D=2.7489 / day), open_bottom=false))

    path = tempname() * ".json"
    @test Agate.Models.export_manifest(path, darwin_manifest) == path
    @test_throws ArgumentError Agate.Models.export_manifest(path, (; schema="x"))
    replayed_from_path = Agate.Models.construct_from_manifest(path; grid=dummy_grid(Float32))
    @test typeof(replayed_from_path) == typeof(darwin_bgc)

    roundtrip_matrix = Float32[0.8 0.2; 0.3 0.7]
    roundtrip_sinking = (P1=0.125f0 / day, D=1.5f0 / day)
    roundtrip_bgc, roundtrip_manifest = Agate.Models.NiPiZD.construct_with_manifest(
        ;
        grid=BoxModelGrid(),
        scalar_type=Float32,
        palatability_matrix=roundtrip_matrix,
        sinking_tracers=roundtrip_sinking,
        open_bottom=false,
    )
    roundtrip_path = tempname() * ".json"
    @test Agate.Models.export_manifest(roundtrip_path, roundtrip_manifest) == roundtrip_path

    roundtrip_json = JSON.parsefile(roundtrip_path)
    @test roundtrip_json["recipe"]["kwargs"]["scalar_type"] == "Float32"
    @test roundtrip_json["recipe"]["kwargs"]["open_bottom"] == false
    sinking_entries = roundtrip_json["recipe"]["kwargs"]["sinking_tracers"]
    @test getindex.(sinking_entries, "name") == ["P1", "D"]
    @test getindex.(sinking_entries, "value") ≈ [0.125f0 / day, 1.5f0 / day]

    matrix_rows = roundtrip_json["recipe"]["kwargs"]["parameters"]["palatability_matrix"]
    decoded_matrix = [matrix_rows[i][j] for i in eachindex(matrix_rows), j in eachindex(matrix_rows[1])]
    @test decoded_matrix ≈ roundtrip_matrix

    replayed_roundtrip = Agate.Models.construct_from_manifest(
        roundtrip_path; grid=BoxModelGrid()
    )
    @test typeof(replayed_roundtrip) == typeof(roundtrip_bgc)
    @test required_biogeochemical_tracers(replayed_roundtrip) ==
        required_biogeochemical_tracers(roundtrip_bgc)
    @test !isnothing(replayed_roundtrip.sinking_velocities)
    @test hasproperty(replayed_roundtrip.sinking_velocities, :P1)
    @test hasproperty(replayed_roundtrip.sinking_velocities, :D)
    @test replayed_roundtrip.parameters.palatability_matrix ≈ roundtrip_matrix

    @test_throws ArgumentError Agate.Models.construct_from_manifest(Dict{String,Any}())
    @test_throws ArgumentError Agate.Models.construct_from_manifest(Dict{String,Any}("schema" => "agate.construction_manifest.v2"))
    @test_throws ArgumentError Agate.Models.construct_from_manifest(Dict{String,Any}("schema" => "agate.construction_manifest.v1"))

    missing_type_manifest = deepcopy(darwin_manifest)
    delete!(missing_type_manifest["recipe"], "type")
    @test_throws ArgumentError Agate.Models.construct_from_manifest(missing_type_manifest; grid=dummy_grid(Float32))

    unsupported_type_manifest = deepcopy(darwin_manifest)
    unsupported_type_manifest["recipe"]["type"] = "factory_graph"
    @test_throws ArgumentError Agate.Models.construct_from_manifest(unsupported_type_manifest; grid=dummy_grid(Float32))

    missing_family_manifest = deepcopy(darwin_manifest)
    delete!(missing_family_manifest["recipe"], "family")
    @test_throws ArgumentError Agate.Models.construct_from_manifest(
        missing_family_manifest; grid=dummy_grid(Float32)
    )

    unsupported_family_manifest = deepcopy(darwin_manifest)
    unsupported_family_manifest["recipe"]["family"] = "NPZD"
    @test_throws ArgumentError Agate.Models.construct_from_manifest(
        unsupported_family_manifest; grid=dummy_grid(Float32)
    )

    constructor_only_manifest = deepcopy(darwin_manifest)
    delete!(constructor_only_manifest["recipe"], "family")
    constructor_only_manifest["recipe"]["constructor"] = "Agate.Models.DARWIN.construct"
    @test_throws ArgumentError Agate.Models.construct_from_manifest(
        constructor_only_manifest; grid=dummy_grid(Float32)
    )
end
