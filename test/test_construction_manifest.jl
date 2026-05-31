using JSON
using Oceananigans.Units: day
using OceanBioME: BoxModelGrid
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers

struct UnsupportedManifestValue end

function darwin_manifest_fixture()
    kwargs = (;
        grid=dummy_grid(Float32),
        phyto_size_structure=(n=3, min_esd=1.5, max_esd=20.0, splitting=:log_splitting),
        zoo_size_structure=(n=2, min_esd=20.0, max_esd=100.0, splitting=:log_splitting),
    )

    bgc = Agate.Models.DARWIN.construct(; kwargs...)
    traced_bgc, manifest = Agate.Models.DARWIN.construct_with_manifest(; kwargs...)
    return bgc, traced_bgc, manifest
end

function nipizd_manifest_fixture(; kwargs...)
    bgc = Agate.Models.NiPiZD.construct(; kwargs...)
    traced_bgc, manifest = Agate.Models.NiPiZD.construct_with_manifest(; kwargs...)
    return bgc, traced_bgc, manifest
end

function test_manifest_roundtrip(manifest, expected_bgc; grid)
    replayed = Agate.Manifests.construct_from_manifest(manifest; grid)
    @test typeof(replayed) == typeof(expected_bgc)
    @test required_biogeochemical_tracers(replayed) == required_biogeochemical_tracers(expected_bgc)
    return replayed
end

@testset "Manifest value serialization" begin
    value = Agate.Manifests.Serialization.manifest_value

    @test value(:foo) == "foo"
    @test value((1, :x, NaN)) == Any[1, "x", "NaN"]
    @test value((a=1, b=:x)) == Dict{String,Any}("a" => 1, "b" => "x")
    @test value(Dict(:a => [1, Inf])) == Dict{String,Any}("a" => Any[1, "Inf"])
    @test_throws ArgumentError value(UnsupportedManifestValue())
end

@testset "DARWIN construction manifests" begin
    darwin_bgc, traced_darwin, manifest = darwin_manifest_fixture()

    @test typeof(traced_darwin) == typeof(darwin_bgc)
    @test manifest["model"]["id"] == "DARWIN/default"
    @test manifest["model"]["family"] == "DARWIN"
    @test !haskey(manifest, "replayable")
    @test manifest["resolved"]["scalar_type"] == "Float32"
    @test manifest["resolved"]["sinking"] == Dict{String,Any}(
        "enabled" => false, "tracers" => nothing, "open_bottom" => true
    )
    @test "P3" in manifest["resolved"]["tracers"]
    @test manifest["recipe"]["type"] == "model_constructor"
    @test manifest["recipe"]["family"] == "DARWIN"
    @test !haskey(manifest["recipe"], "constructor")
    @test manifest["recipe"]["kwargs"]["phyto_size_structure"] ==
        manifest["resolved"]["plankton_diameters_by_group"]["P"]
    @test manifest["recipe"]["kwargs"]["zoo_size_structure"] ==
        manifest["resolved"]["plankton_diameters_by_group"]["Z"]
    @test manifest["recipe"]["kwargs"]["parameters"] isa Dict
    @test manifest["recipe"]["kwargs"]["sinking_tracers"] === nothing
    @test manifest["recipe"]["kwargs"]["open_bottom"] == true
    @test manifest["recipe"]["kwargs"]["scalar_type"] == "Float32"
end

@testset "NiPiZD construction manifests" begin
    nipizd_bgc, traced_nipizd, manifest = nipizd_manifest_fixture(; grid=dummy_grid(Float32))

    @test typeof(traced_nipizd) == typeof(nipizd_bgc)
    @test manifest["model"]["id"] == "NiPiZD/default"
    @test manifest["model"]["family"] == "NiPiZD"
    @test !haskey(manifest, "replayable")
    @test manifest["resolved"]["scalar_type"] == "Float32"
    @test haskey(manifest["resolved"], "parameters")
    @test manifest["recipe"]["type"] == "model_constructor"
    @test manifest["recipe"]["family"] == "NiPiZD"
    @test !haskey(manifest["recipe"], "constructor")
    @test haskey(manifest["recipe"]["kwargs"], "parameters")
end

@testset "Sinking manifest entries" begin
    _, manifest = Agate.Models.NiPiZD.construct_with_manifest(
        ;
        grid=BoxModelGrid(),
        sinking_tracers=(P1=0.2551 / day, D=2.7489 / day),
        open_bottom=false,
    )

    @test manifest["resolved"]["has_sinking_velocities"]
    @test manifest["resolved"]["sinking"]["enabled"]
    @test manifest["resolved"]["sinking"]["open_bottom"] == false
    @test manifest["resolved"]["sinking"]["tracers"] == Dict{String,Any}(
        "P1" => 0.2551 / day, "D" => 2.7489 / day
    )
    @test manifest["recipe"]["kwargs"]["sinking_tracers"] == Any[
        Dict{String,Any}("name" => "P1", "value" => 0.2551 / day),
        Dict{String,Any}("name" => "D", "value" => 2.7489 / day),
    ]
    @test manifest["recipe"]["kwargs"]["open_bottom"] == false
end

@testset "Manifest replay" begin
    darwin_bgc, _, darwin_manifest = darwin_manifest_fixture()
    nipizd_bgc, _, nipizd_manifest = nipizd_manifest_fixture(; grid=dummy_grid(Float32))

    test_manifest_roundtrip(darwin_manifest, darwin_bgc; grid=dummy_grid(Float32))
    test_manifest_roundtrip(nipizd_manifest, nipizd_bgc; grid=dummy_grid(Float32))

    sinking_bgc, sinking_manifest = Agate.Models.NiPiZD.construct_with_manifest(
        ; grid=BoxModelGrid(), sinking_tracers=(P1=0.2551 / day, D=2.7489 / day), open_bottom=false
    )
    replayed_sinking = Agate.Manifests.construct_from_manifest(sinking_manifest; grid=BoxModelGrid())
    @test typeof(replayed_sinking) == typeof(sinking_bgc)
end

@testset "Manifest file IO" begin
    darwin_bgc, _, manifest = darwin_manifest_fixture()
    path = tempname() * ".json"

    @test Agate.Manifests.export_manifest(path, manifest) == path
    @test_throws ArgumentError Agate.Manifests.export_manifest(path, (; schema="x"))

    replayed = Agate.Manifests.construct_from_manifest(path; grid=dummy_grid(Float32))
    @test typeof(replayed) == typeof(darwin_bgc)
end

@testset "Manifest roundtrip values" begin
    matrix = Float32[0.8 0.2; 0.3 0.7]
    sinking = (P1=0.125f0 / day, D=1.5f0 / day)
    bgc, manifest = Agate.Models.NiPiZD.construct_with_manifest(
        ;
        grid=BoxModelGrid(),
        scalar_type=Float32,
        palatability_matrix=matrix,
        sinking_tracers=sinking,
        open_bottom=false,
    )

    path = tempname() * ".json"
    @test Agate.Manifests.export_manifest(path, manifest) == path

    json = JSON.parsefile(path)
    @test json["recipe"]["kwargs"]["scalar_type"] == "Float32"
    @test json["recipe"]["kwargs"]["open_bottom"] == false

    sinking_entries = json["recipe"]["kwargs"]["sinking_tracers"]
    @test getindex.(sinking_entries, "name") == ["P1", "D"]
    @test getindex.(sinking_entries, "value") ≈ [0.125f0 / day, 1.5f0 / day]

    matrix_rows = json["recipe"]["kwargs"]["parameters"]["palatability_matrix"]
    decoded_matrix = [matrix_rows[i][j] for i in eachindex(matrix_rows), j in eachindex(matrix_rows[1])]
    @test decoded_matrix ≈ matrix

    replayed = Agate.Manifests.construct_from_manifest(path; grid=BoxModelGrid())
    @test typeof(replayed) == typeof(bgc)
    @test required_biogeochemical_tracers(replayed) == required_biogeochemical_tracers(bgc)
    @test !isnothing(replayed.sinking_velocities)
    @test hasproperty(replayed.sinking_velocities, :P1)
    @test hasproperty(replayed.sinking_velocities, :D)
    @test replayed.parameters.palatability_matrix ≈ matrix
end

@testset "Manifest validation errors" begin
    _, _, manifest = darwin_manifest_fixture()

    @test_throws ArgumentError Agate.Manifests.construct_from_manifest(Dict{String,Any}())
    @test_throws ArgumentError Agate.Manifests.construct_from_manifest(
        Dict{String,Any}("schema" => "agate.construction_manifest.v2")
    )
    @test_throws ArgumentError Agate.Manifests.construct_from_manifest(
        Dict{String,Any}("schema" => "agate.construction_manifest.v1")
    )

    missing_type_manifest = deepcopy(manifest)
    delete!(missing_type_manifest["recipe"], "type")
    @test_throws ArgumentError Agate.Manifests.construct_from_manifest(
        missing_type_manifest; grid=dummy_grid(Float32)
    )

    unsupported_type_manifest = deepcopy(manifest)
    unsupported_type_manifest["recipe"]["type"] = "factory_graph"
    @test_throws ArgumentError Agate.Manifests.construct_from_manifest(
        unsupported_type_manifest; grid=dummy_grid(Float32)
    )

    missing_family_manifest = deepcopy(manifest)
    delete!(missing_family_manifest["recipe"], "family")
    @test_throws ArgumentError Agate.Manifests.construct_from_manifest(
        missing_family_manifest; grid=dummy_grid(Float32)
    )

    unsupported_family_manifest = deepcopy(manifest)
    unsupported_family_manifest["recipe"]["family"] = "NPZD"
    @test_throws ArgumentError Agate.Manifests.construct_from_manifest(
        unsupported_family_manifest; grid=dummy_grid(Float32)
    )

    constructor_only_manifest = deepcopy(manifest)
    delete!(constructor_only_manifest["recipe"], "family")
    constructor_only_manifest["recipe"]["constructor"] = "Agate.Models.DARWIN.construct"
    @test_throws ArgumentError Agate.Manifests.construct_from_manifest(
        constructor_only_manifest; grid=dummy_grid(Float32)
    )
end
