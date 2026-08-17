using Oceananigans.Units: day
using OceanBioME: BoxModelGrid
using Oceananigans.Biogeochemistry: biogeochemical_drift_velocity, required_biogeochemical_tracers
using Agate.Manifests: construct_from_manifest, export_manifest

function Agate.Manifests.construct_from_manifest(
    ::Val{:TestModelSetup}, setup::AbstractDict; grid=nothing, arch=nothing
)
    return (; value=setup["kwargs"]["value"], grid, arch)
end

function test_reconstructed_model(reconstructed, expected)
    @test typeof(reconstructed) == typeof(expected)
    @test required_biogeochemical_tracers(reconstructed) == required_biogeochemical_tracers(expected)
end

function test_exact_parameters(actual, expected)
    @test typeof(actual) == typeof(expected)
    for key in keys(expected)
        actual_value = getproperty(actual, key)
        expected_value = getproperty(expected, key)
        if key === :interactions
            @test actual_value == expected_value
        else
            @test isequal(actual_value, expected_value)
        end
    end
end

function test_manifest_argument_error(edit, setup, message)
    invalid = deepcopy(setup)
    edit(invalid)
    err = try
        construct_from_manifest(invalid)
        nothing
    catch err
        err
    end
    @test err isa ArgumentError
    err isa ArgumentError && @test occursin(message, sprint(showerror, err))
end

@testset "NiPiZD Float32 model setup round trip" begin
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

    @test setup["agate"]["version"] == string(Base.pkgversion(Agate))

    path = tempname() * ".json"
    @test export_manifest(path, setup) == path

    reconstructed = construct_from_manifest(path; grid)
    test_reconstructed_model(reconstructed, bgc)

    test_exact_parameters(reconstructed.parameters, bgc.parameters)
    @test !isnothing(reconstructed.sinking_velocities)

    reconstructed_P_1 = biogeochemical_drift_velocity(reconstructed, Val(:P_1)).w.data[1, 1, 1]
    reconstructed_D = biogeochemical_drift_velocity(reconstructed, Val(:D)).w.data[1, 1, 1]

    @test reconstructed_P_1 == biogeochemical_drift_velocity(bgc, Val(:P_1)).w.data[1, 1, 1]
    @test reconstructed_D == biogeochemical_drift_velocity(bgc, Val(:D)).w.data[1, 1, 1]
    @test reconstructed_P_1 == -sinking_tracers.P_1
    @test reconstructed_D == -sinking_tracers.D
end

@testset "NiPiZD non-finite model setup round trip" begin
    bgc, setup = Agate.Models.NiPiZD.construct_with_manifest(
        ;
        scalar_type=Float32,
        parameters=(;
            detritus_remineralization=Float32(NaN),
            linear_mortality=Float32[Inf, -Inf, 0, 0],
        ),
        palatability_matrix=Float32[NaN Inf; -Inf 0.5],
    )

    path = tempname() * ".json"
    export_manifest(path, setup)
    reconstructed = construct_from_manifest(path)
    test_reconstructed_model(reconstructed, bgc)

    test_exact_parameters(reconstructed.parameters, bgc.parameters)
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

    test_exact_parameters(reconstructed.parameters, bgc.parameters)
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
    test_exact_parameters(reconstructed.parameters, bgc.parameters)
    @test Agate.Introspection.plankton_diameters(reconstructed) ==
          Agate.Introspection.plankton_diameters(bgc)
end

@testset "Model setup family dispatch" begin
    setup = Dict{String,Any}(
        "schema" => "agate.model_setup.v1",
        "model" => Dict{String,Any}("family" => "TestModelSetup"),
        "kwargs" => Dict{String,Any}("value" => 3),
    )

    @test construct_from_manifest(setup; grid=:grid, arch=:arch) ==
          (value=3, grid=:grid, arch=:arch)

    setup["model"]["family"] = "UnknownModel"
    @test_throws ArgumentError construct_from_manifest(setup)
end

@testset "Model setup serialization validation" begin
    @test_throws ArgumentError Agate.Manifests.Serialization.manifest_value(:unsupported)
end

@testset "Model setup reader validation" begin
    _, setup = Agate.Models.NiPiZD.construct_with_manifest(;
        palatability_matrix=[0.8 0.2; 0.3 0.7]
    )

    test_manifest_argument_error(setup, "Model setup has unsupported field") do invalid
        invalid["unexpected"] = true
    end
    test_manifest_argument_error(setup, "agate.model_setup.v2") do invalid
        invalid["schema"] = "agate.model_setup.v2"
    end
    test_manifest_argument_error(setup, "Model setup kwargs has unsupported field") do invalid
        invalid["kwargs"]["unexpected"] = true
    end
    test_manifest_argument_error(setup, "missing required field \"scalar_type\"") do invalid
        delete!(invalid["kwargs"], "scalar_type")
    end
    test_manifest_argument_error(setup, "parameters must be an object") do invalid
        invalid["kwargs"]["parameters"] = Any[]
    end
    test_manifest_argument_error(setup, "size_structure.phytoplankton[1]") do invalid
        invalid["kwargs"]["size_structure"]["phytoplankton"][1]["unexpected"] = true
    end
    test_manifest_argument_error(setup, "sinking_tracers[1]") do invalid
        invalid["kwargs"]["sinking_tracers"] = [Dict("name" => "D")]
    end
    test_manifest_argument_error(setup, "ambiguous empty array") do invalid
        invalid["kwargs"]["parameters"]["palatability_matrix"] = Any[]
    end
    test_manifest_argument_error(setup, "matrix must be rectangular") do invalid
        invalid["kwargs"]["parameters"]["palatability_matrix"] = [[1.0, 2.0], [3.0]]
    end
end
