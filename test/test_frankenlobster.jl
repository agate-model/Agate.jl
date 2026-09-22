using Test
using Adapt
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers

const FrankenLOBSTER = Agate.Models.FrankenLOBSTER

@testset "FrankenLOBSTER construction boundary" begin
    plankton = FrankenLOBSTER._construct_plankton(; grid=dummy_grid(Float32))
    ownership = (
        required_biogeochemical_tracers(plankton),
        FrankenLOBSTER.external_tracers(plankton),
    )

    @test ownership == ((:P_1, :P_2, :Z_1, :Z_2, :B_1), (:NO₃, :NH₄, :DOM))
    @test required_biogeochemical_tracers(plankton.runtime) ==
        (:NO₃, :NH₄, :DOM, :P_1, :P_2, :Z_1, :Z_2, :B_1)
    @test plankton.runtime.metadata.plankton_diameters ==
        (0.6f0, 1.2f0, 6.0f0, 12.0f0, 0.6f0)

    adapted = Adapt.adapt(identity, plankton)
    @test (
        required_biogeochemical_tracers(adapted),
        FrankenLOBSTER.external_tracers(adapted),
    ) == ownership
end

@testset "FrankenLOBSTER arbitrary P/Z/B realization" begin
    plankton = FrankenLOBSTER._construct_plankton(;
        size_structure=(
            phytoplankton=(pico=[0.5], nano=[2.0]),
            zooplankton=(micro=[8.0], meso=[20.0]),
            bacterioplankton=(heterotroph=[0.8, 1.6],),
        ),
        grid=dummy_grid(Float64),
    )

    @test required_biogeochemical_tracers(plankton) ==
        (:nano_1, :pico_1, :meso_1, :micro_1, :heterotroph_1, :heterotroph_2)
end
