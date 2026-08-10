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
    sinking_tracers = (P_1=0.125f0 / day, D=1.5f0 / day)

    bgc, setup = Agate.Models.NiPiZD.construct_with_manifest(
        ;
        grid,
        scalar_type=Float32,
        palatability_matrix,
        sinking_tracers,
    )

    path = tempname() * ".json"
    @test export_manifest(path, setup) == path

    reconstructed = construct_from_manifest(path; grid)
    test_reconstructed_model(reconstructed, bgc)

    @test reconstructed.parameters.palatability_matrix ≈ palatability_matrix
    @test !isnothing(reconstructed.sinking_velocities)

    reconstructed_P_1 = biogeochemical_drift_velocity(reconstructed, Val(:P_1)).w.data[1, 1, 1]
    reconstructed_D = biogeochemical_drift_velocity(reconstructed, Val(:D)).w.data[1, 1, 1]

    @test reconstructed_P_1 ≈ biogeochemical_drift_velocity(bgc, Val(:P_1)).w.data[1, 1, 1]
    @test reconstructed_D ≈ biogeochemical_drift_velocity(bgc, Val(:D)).w.data[1, 1, 1]
    @test reconstructed_P_1 ≈ -sinking_tracers.P_1
    @test reconstructed_D ≈ -sinking_tracers.D
end

@testset "NiPiZD named model setup round trip" begin
    size_structure = (;
        phytoplankton=(diat=[2.0, 5.0], dino=[10.0]),
        zooplankton=(microzoo=[30.0, 60.0], mesozoo=[100.0]),
    )
    bgc, setup = Agate.Models.NiPiZD.construct_with_manifest(;
        size_structure,
        parameters=(;
            maximum_growth_rate=(diat_2=1.25 / day,),
            specificity=(microzoo_1=2.0,),
        ),
    )

    path = tempname() * ".json"
    export_manifest(path, setup)
    reconstructed = construct_from_manifest(path)
    test_reconstructed_model(reconstructed, bgc)

    @test Agate.Introspection.plankton_groups(reconstructed) ==
          Agate.Introspection.plankton_groups(bgc)
    @test Agate.Introspection.plankton_diameters(reconstructed) ==
          Agate.Introspection.plankton_diameters(bgc)

    active(model) = Agate.Runtime.active_parameters(model;
        maximum_growth_rate=(:diat_2,),
        interactions=(; palatability=((:microzoo_1, :diat_1),)),
    )
    @test active(reconstructed).labels == active(bgc).labels
    @test active(reconstructed).values == active(bgc).values
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
