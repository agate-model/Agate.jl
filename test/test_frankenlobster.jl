using Test
using Adapt
using Oceananigans.Biogeochemistry:
    required_biogeochemical_auxiliary_fields,
    required_biogeochemical_tracers

using OceanBioME.Models.NutrientsPlanktonDetritusModels:
    InstantRemineralisationDetritus,
    NutrientsPlanktonDetritus
using OceanBioME.Models.NutrientsPlanktonDetritusModels.NutrientsModels:
    Nutrients,
    NitrateAmmonia
using OceanBioME.Models.NutrientsPlanktonDetritusModels:
    nutrient_uptake,
    solid_waste

const FrankenLOBSTER = Agate.Models.FrankenLOBSTER

_cell(value) = fill(value, 1, 1, 1)

function _frankenlobster_fields(; NO₃=1.0, NH₄=1.0, DOM=0.0,
                                P_1=0.0, P_2=0.0, Z_1=0.0, Z_2=0.0, B_1=0.0)
    return (
        NO₃=_cell(NO₃),
        NH₄=_cell(NH₄),
        DOM=_cell(DOM),
        P_1=_cell(P_1),
        P_2=_cell(P_2),
        Z_1=_cell(Z_1),
        Z_2=_cell(Z_2),
        B_1=_cell(B_1),
    )
end

function _frankenlobster_npd(plankton; nitrification_rate=0.0)
    nutrients = Nutrients(; nitrogen=NitrateAmmonia(; nitrification_rate))
    return NutrientsPlanktonDetritus{Float64}(
        nutrients, plankton, InstantRemineralisationDetritus(), nothing, nothing
    )
end

@testset "FrankenLOBSTER construction boundary" begin
    plankton = FrankenLOBSTER._construct_plankton(; grid=dummy_grid(Float32))
    ownership = (
        required_biogeochemical_tracers(plankton),
        FrankenLOBSTER.external_tracers(plankton),
        FrankenLOBSTER.exchange_tracers(plankton),
    )

    @test ownership == (
        (:P_1, :P_2, :Z_1, :Z_2, :B_1),
        (:NO₃, :NH₄, :DOM),
        (:solid_waste,),
    )
    @test required_biogeochemical_tracers(plankton.runtime) ==
        (:NO₃, :NH₄, :DOM, :solid_waste, :P_1, :P_2, :Z_1, :Z_2, :B_1)
    @test required_biogeochemical_auxiliary_fields(plankton) == (:PAR,)
    @test plankton.runtime.metadata.plankton_diameters ==
        (0.6f0, 1.2f0, 6.0f0, 12.0f0, 0.6f0)

    adapted = Adapt.adapt(identity, plankton)
    @test (
        required_biogeochemical_tracers(adapted),
        FrankenLOBSTER.external_tracers(adapted),
        FrankenLOBSTER.exchange_tracers(adapted),
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

@testset "FrankenLOBSTER NPD nitrate/ammonium growth bridge" begin
    plankton = FrankenLOBSTER._construct_plankton(;
        grid=dummy_grid(Float64),
        parameters=(
            maximum_growth_rate=(P_1=1.0, P_2=1.0),
            nitrate_half_saturation=(P_1=1.0, P_2=1.0),
            ammonium_half_saturation=(P_1=1.0, P_2=1.0),
            light_half_saturation=(P_1=1.0, P_2=1.0),
            nitrate_ammonia_inhibition=log(2.0),
            phytoplankton_mortality_rate=(P_1=0.0, P_2=0.0),
            zooplankton_mortality_rate=(Z_1=0.0, Z_2=0.0),
            maximum_predation_rate=(Z_1=0.0, Z_2=0.0),
        ),
    )
    bgc = _frankenlobster_npd(plankton; nitrification_rate=0.1)
    fields = _frankenlobster_fields(; P_1=2.0)
    auxiliary_fields = (PAR=_cell(1.0),)
    clock = (; time=0.0)
    grid = dummy_grid(Float64)

    @test bgc(1, 1, 1, grid, Val(:P_1), clock, fields, auxiliary_fields) ≈ 0.75
    @test nutrient_uptake(
        1, 1, 1, grid, Val(:NO₃), plankton, bgc, fields, auxiliary_fields
    ) ≈ 0.25
    @test nutrient_uptake(
        1, 1, 1, grid, Val(:NH₄), plankton, bgc, fields, auxiliary_fields
    ) ≈ 0.5
    @test bgc(1, 1, 1, grid, Val(:NO₃), clock, fields, auxiliary_fields) ≈ -0.15
    @test bgc(1, 1, 1, grid, Val(:NH₄), clock, fields, auxiliary_fields) ≈ -0.6
end

@testset "FrankenLOBSTER P/Z losses route to NPD solid waste" begin
    plankton = FrankenLOBSTER._construct_plankton(;
        grid=dummy_grid(Float64),
        parameters=(
            maximum_growth_rate=(P_1=0.0, P_2=0.0),
            phytoplankton_mortality_rate=(P_1=0.25, P_2=0.0),
            zooplankton_mortality_rate=(Z_1=0.0, Z_2=0.0),
            maximum_predation_rate=(Z_1=0.0, Z_2=0.0),
        ),
    )
    bgc = _frankenlobster_npd(plankton)
    fields = _frankenlobster_fields(; P_1=2.0)
    auxiliary_fields = (PAR=_cell(1.0),)
    clock = (; time=0.0)
    grid = dummy_grid(Float64)

    p_loss = bgc(1, 1, 1, grid, Val(:P_1), clock, fields, auxiliary_fields)
    waste = solid_waste(1, 1, 1, grid, plankton, bgc, fields, auxiliary_fields)
    @test p_loss ≈ -1.0
    @test waste ≈ 1.0
    @test p_loss + waste ≈ 0.0 atol=1e-14
end

@testset "FrankenLOBSTER Z grazing conserves living transfer and waste" begin
    plankton = FrankenLOBSTER._construct_plankton(;
        grid=dummy_grid(Float64),
        parameters=(
            maximum_growth_rate=(P_1=0.0, P_2=0.0),
            phytoplankton_mortality_rate=(P_1=0.0, P_2=0.0),
            zooplankton_mortality_rate=(Z_1=0.0, Z_2=0.0),
            maximum_predation_rate=(Z_1=1.0, Z_2=0.0),
            grazing_half_saturation=(Z_1=1.0, Z_2=1.0),
            palatability_matrix=[1.0 0.0; 0.0 0.0],
            assimilation_matrix=[0.5 0.0; 0.0 0.0],
        ),
    )
    bgc = _frankenlobster_npd(plankton)
    fields = _frankenlobster_fields(; P_1=2.0, Z_1=1.0)
    auxiliary_fields = (PAR=_cell(1.0),)
    clock = (; time=0.0)
    grid = dummy_grid(Float64)

    p = bgc(1, 1, 1, grid, Val(:P_1), clock, fields, auxiliary_fields)
    z = bgc(1, 1, 1, grid, Val(:Z_1), clock, fields, auxiliary_fields)
    waste = solid_waste(1, 1, 1, grid, plankton, bgc, fields, auxiliary_fields)

    @test p ≈ -2 / 3
    @test z ≈ 1 / 3
    @test waste ≈ 1 / 3
    @test p + z + waste ≈ 0.0 atol=1e-14
end
