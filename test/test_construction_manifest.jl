using Adapt
using JSON
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
    @test "P3" in darwin_manifest["resolved"]["tracers"]

    nipizd_bgc = Agate.Models.NiPiZD.construct(; grid=dummy_grid(Float32))
    traced_nipizd, nipizd_manifest = Agate.Models.NiPiZD.construct_with_manifest(
        ; grid=dummy_grid(Float32)
    )

    @test typeof(traced_nipizd) == typeof(nipizd_bgc)
    @test nipizd_manifest["model"]["id"] == "NiPiZD/default"
    @test nipizd_manifest["model"]["family"] == "NiPiZD"
    @test nipizd_manifest["resolved"]["scalar_type"] == "Float32"
    @test haskey(nipizd_manifest["resolved"], "parameters")
end
