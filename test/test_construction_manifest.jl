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

@testset "Construction manifest" begin
    id = Agate.Models.ModelId(:DARWIN, :citation2026, :A)
    spec = Agate.Models.variant(
        id;
        phyto_size_structure=(n=3, min_esd=2, max_esd=10, splitting=:log_splitting),
        zoo_size_structure=(n=2, min_esd=20, max_esd=100, splitting=:linear_splitting),
    )

    bgc = Agate.Models.construct(spec; grid=dummy_grid(Float32))
    @test !hasproperty(bgc, :agate_metadata)

    traced_bgc, manifest = Agate.Models.construct_with_manifest(spec; grid=dummy_grid(Float32))

    @test typeof(traced_bgc) == typeof(bgc)
    @test manifest["schema"] == "agate.construction_manifest.v1"
    @test manifest["model"]["id"] == "DARWIN/citation2026/A"
    @test !haskey(manifest, "requested")
    @test haskey(manifest["agate"], "version")
    @test !isempty(manifest["resolved"]["tracers"])
    @test manifest["resolved"]["scalar_type"] == "Float32"
    @test manifest["resolved"]["auxiliary_fields"] == Any["PAR"]
    @test haskey(manifest["resolved"], "parameters")
    @test haskey(manifest["resolved"]["parameters"], "palatability_matrix")
    @test manifest["resolved"]["parameters"]["palatability_matrix"]["shape"] == "matrix"
    @test manifest["resolved"]["parameters"]["palatability_matrix"]["value"] isa Vector

    path = tempname() * ".json"
    @test Agate.Models.export_manifest(path, manifest) == path
    @test_throws ArgumentError Agate.Models.export_manifest(path, (; schema="x"))
    @test isfile(path)
    parsed = JSON.parsefile(path)
    @test parsed["schema"] == "agate.construction_manifest.v1"
    @test parsed["model"]["id"] == "DARWIN/citation2026/A"

    @test !isempty(required_biogeochemical_tracers(traced_bgc))
    adapted = Adapt.adapt(identity, traced_bgc)
    @test typeof(adapted) == typeof(traced_bgc)
    @test !hasproperty(adapted, :agate_metadata)
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
    @test darwin_manifest["resolved"]["scalar_type"] == "Float32"
    @test darwin_manifest["resolved"]["sinking"] == Dict{String,Any}(
        "enabled" => false, "tracers" => nothing, "open_bottom" => true
    )
    @test "P3" in darwin_manifest["resolved"]["tracers"]
    @test darwin_manifest["recipe"]["constructor"] == "Agate.Models.DARWIN.construct"
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
    @test nipizd_manifest["resolved"]["scalar_type"] == "Float32"
    @test haskey(nipizd_manifest["resolved"], "parameters")
    @test nipizd_manifest["recipe"]["constructor"] == "Agate.Models.NiPiZD.construct"
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
    Agate.Models.export_manifest(path, darwin_manifest)
    replayed_from_path = Agate.Models.construct_from_manifest(path; grid=dummy_grid(Float32))
    @test typeof(replayed_from_path) == typeof(darwin_bgc)

    @test_throws ArgumentError Agate.Models.construct_from_manifest(Dict{String,Any}("schema" => "agate.construction_manifest.v1"))
end
