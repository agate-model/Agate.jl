using Oceananigans.Units: day
using OceanBioME: BoxModelGrid
using Oceananigans.Biogeochemistry: biogeochemical_drift_velocity, required_biogeochemical_tracers
using Agate.Manifests: construct_from_manifest, export_manifest

function test_reconstructed_model(reconstructed, expected)
    @test typeof(reconstructed) == typeof(expected)
    @test required_biogeochemical_tracers(reconstructed) == required_biogeochemical_tracers(expected)
end

@testset "NiPiZD model setup export and import" begin
    grid = BoxModelGrid()
    palatability_matrix = Float32[0.8 0.2; 0.3 0.7]
    sinking_tracers = (P1=0.125f0 / day, D=1.5f0 / day)

    bgc, setup = Agate.Models.NiPiZD.construct_with_manifest(
        ;
        grid,
        scalar_type=Float32,
        palatability_matrix,
        sinking_tracers,
        open_bottom=false,
    )

    path = tempname() * ".json"
    @test export_manifest(path, setup) == path

    reconstructed = construct_from_manifest(path; grid)
    test_reconstructed_model(reconstructed, bgc)

    @test reconstructed.parameters.palatability_matrix ≈ palatability_matrix
    @test !isnothing(reconstructed.sinking_velocities)
    @test biogeochemical_drift_velocity(reconstructed, Val(:P1)).w.data[1, 1, 1] ≈ -sinking_tracers.P1
    @test biogeochemical_drift_velocity(reconstructed, Val(:D)).w.data[1, 1, 1] ≈ -sinking_tracers.D
end

@testset "DARWIN model setup import" begin
    grid = dummy_grid(Float32)

    bgc, setup = Agate.Models.DARWIN.construct_with_manifest(
        ;
        grid,
        phyto_size_structure=(n=3, min_esd=1.5, max_esd=20.0, splitting=:log_splitting),
        zoo_size_structure=(n=2, min_esd=20.0, max_esd=100.0, splitting=:log_splitting),
    )

    reconstructed = construct_from_manifest(setup; grid)
    test_reconstructed_model(reconstructed, bgc)
end
