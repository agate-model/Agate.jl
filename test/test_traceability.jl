using Adapt
using JSON
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers

@testset "Construction traceability" begin
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
    @test haskey(manifest["resolved"], "parameters")
    @test haskey(manifest["resolved"]["parameters"], "palatability_matrix")
    @test manifest["resolved"]["parameters"]["palatability_matrix"]["shape"] == "matrix"
    @test manifest["resolved"]["parameters"]["palatability_matrix"]["value"] isa Vector

    path = tempname() * ".json"
    @test Agate.Traceability.export_manifest(path, manifest) == path
    @test isfile(path)
    parsed = JSON.parsefile(path)
    @test parsed["schema"] == "agate.construction_manifest.v1"
    @test parsed["model"]["id"] == "DARWIN/citation2026/A"

    @test !isempty(required_biogeochemical_tracers(traced_bgc))
    adapted = Adapt.adapt(identity, traced_bgc)
    @test typeof(adapted) == typeof(traced_bgc)
    @test !hasproperty(adapted, :agate_metadata)
end
