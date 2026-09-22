using Test
using Adapt
using Oceananigans.Architectures: CPU
using Oceananigans.Grids: RectilinearGrid
using Oceananigans.Biogeochemistry:
    required_biogeochemical_auxiliary_fields,
    required_biogeochemical_tracers

using OceanBioME.Models.NutrientsPlanktonDetritusModels:
    DissolvedParticulate,
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

function _frankenlobster_fields(; NO₃=1.0, NH₄=1.0, DOM=0.0, sPOM=0.0, bPOM=0.0,
                                P_1=0.0, P_2=0.0, Z_1=0.0, Z_2=0.0, B_1=0.0)
    return (
        NO₃=_cell(NO₃),
        NH₄=_cell(NH₄),
        DOM=_cell(DOM),
        sPOM=_cell(sPOM),
        bPOM=_cell(bPOM),
        P_1=_cell(P_1),
        P_2=_cell(P_2),
        Z_1=_cell(Z_1),
        Z_2=_cell(Z_2),
        B_1=_cell(B_1),
    )
end

function _frankenlobster_npd(
    plankton; nitrification_rate=0.0, detritus=InstantRemineralisationDetritus()
)
    nutrients = Nutrients(; nitrogen=NitrateAmmonia(; nitrification_rate))
    return NutrientsPlanktonDetritus{Float64}(nutrients, plankton, detritus, nothing, nothing)
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
        (:solid_waste, :inorganic_waste),
    )
    @test required_biogeochemical_tracers(plankton.runtime) ==
        (:NO₃, :NH₄, :DOM, :solid_waste, :inorganic_waste, :P_1, :P_2, :Z_1, :Z_2, :B_1)
    @test required_biogeochemical_auxiliary_fields(plankton) == (:PAR,)
    @test plankton.runtime.metadata.plankton_diameters ==
        (0.6f0, 1.2f0, 6.0f0, 12.0f0, 0.6f0)

    nitrate_K = plankton.runtime.parameters.nitrate_half_saturation
    ammonium_K = plankton.runtime.parameters.ammonium_half_saturation
    @test ammonium_K ≈ 0.5f0 .* nitrate_K
    @test nitrate_K[1] < nitrate_K[2]

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
    @test size(plankton.runtime.parameters.palatability_matrix) == (2, 4)
    @test size(plankton.runtime.parameters.assimilation_matrix) == (2, 4)
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
            palatability_matrix=[1.0 0.0 0.0; 0.0 0.0 0.0],
            assimilation_matrix=[0.5 0.0 0.0; 0.0 0.0 0.0],
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


@testset "FrankenLOBSTER Z shares one ingestion capacity across P and B prey" begin
    plankton = FrankenLOBSTER._construct_plankton(;
        grid=dummy_grid(Float64),
        parameters=(
            maximum_growth_rate=(P_1=0.0, P_2=0.0),
            phytoplankton_mortality_rate=(P_1=0.0, P_2=0.0),
            zooplankton_mortality_rate=(Z_1=0.0, Z_2=0.0),
            bacterioplankton_mortality_rate=(B_1=0.0,),
            bacterial_maximum_uptake_rate=(B_1=0.0,),
            maximum_predation_rate=(Z_1=1.0, Z_2=0.0),
            grazing_half_saturation=(Z_1=1.0, Z_2=1.0),
            palatability_matrix=[1.0 0.0 1.0; 0.0 0.0 0.0],
            assimilation_matrix=[0.5 0.0 0.25; 0.0 0.0 0.0],
        ),
    )
    bgc = _frankenlobster_npd(plankton)
    fields = _frankenlobster_fields(; P_1=1.0, B_1=1.0, Z_1=1.0)
    auxiliary_fields = (PAR=_cell(1.0),)
    clock = (; time=0.0)
    grid = dummy_grid(Float64)

    p = bgc(1, 1, 1, grid, Val(:P_1), clock, fields, auxiliary_fields)
    b = bgc(1, 1, 1, grid, Val(:B_1), clock, fields, auxiliary_fields)
    z = bgc(1, 1, 1, grid, Val(:Z_1), clock, fields, auxiliary_fields)
    waste = solid_waste(1, 1, 1, grid, plankton, bgc, fields, auxiliary_fields)

    @test [p, b, z, waste] ≈ [-1 / 3, -1 / 3, 1 / 4, 5 / 12]
    @test -(p + b) ≈ 2 / 3
    @test p + b + z + waste ≈ 0.0 atol=1e-14
end


@testset "FrankenLOBSTER bacterial size meta-traits" begin
    plankton = FrankenLOBSTER._construct_plankton(;
        size_structure=(
            phytoplankton=(P=[0.6, 1.2],),
            zooplankton=(Z=[6.0, 12.0],),
            bacterioplankton=(B=[0.4, 0.8],),
        ),
        grid=dummy_grid(Float64),
    )

    mu = plankton.runtime.parameters.bacterial_maximum_uptake_rate
    K = plankton.runtime.parameters.bacterial_dom_half_saturation
    volume(d) = pi / 6 * d^3
    expected_mu = [1.836 / 86400 * volume(d)^0.28 for d in (0.4, 0.8)]
    expected_K = [0.04284 * volume(d)^0.65 for d in (0.4, 0.8)]

    @test mu ≈ expected_mu
    @test vec(K) ≈ expected_K
    @test mu[1] < mu[2] && K[1, 1] < K[2, 1]
end

@testset "FrankenLOBSTER DOM uptake closes through bacterial growth and regeneration" begin
    grid = RectilinearGrid(CPU(); size=(1, 1, 1), extent=(1, 1, 1))
    plankton = FrankenLOBSTER._construct_plankton(;
        grid,
        parameters=(
            maximum_growth_rate=(P_1=0.0, P_2=0.0),
            phytoplankton_mortality_rate=(P_1=0.0, P_2=0.0),
            zooplankton_mortality_rate=(Z_1=0.0, Z_2=0.0),
            maximum_predation_rate=(Z_1=0.0, Z_2=0.0),
            bacterial_maximum_uptake_rate=(B_1=2.0,),
            bacterial_dom_half_saturation=reshape([1.0], 1, 1),
            bacterial_substrate_preference=reshape([1.0], 1, 1),
            bacterial_assimilation=reshape([0.25], 1, 1),
            bacterioplankton_mortality_rate=(B_1=0.0,),
        ),
    )
    detritus = DissolvedParticulate(
        grid;
        dissolved_remineralisation_rate=0.0,
        particulate_remineralisation_rate=(0.0, 0.0),
        sinking_speeds=(0.0, 0.0),
    )
    bgc = _frankenlobster_npd(plankton; detritus)
    fields = _frankenlobster_fields(; DOM=3.0, B_1=2.0)
    auxiliary_fields = (PAR=_cell(1.0),)
    clock = (; time=0.0)

    b_gain = bgc(1, 1, 1, grid, Val(:B_1), clock, fields, auxiliary_fields)
    dom = bgc(1, 1, 1, grid, Val(:DOM), clock, fields, auxiliary_fields)
    nh4 = bgc(1, 1, 1, grid, Val(:NH₄), clock, fields, auxiliary_fields)

    @test [dom, b_gain, nh4] ≈ [-3.0, 0.75, 2.25]
    @test dom + b_gain + nh4 ≈ 0.0 atol=1e-14
    @test bgc(1, 1, 1, grid, Val(:sPOM), clock, fields, auxiliary_fields) == 0.0
    @test bgc(1, 1, 1, grid, Val(:bPOM), clock, fields, auxiliary_fields) == 0.0

    zero_fields = _frankenlobster_fields(; DOM=0.0, B_1=2.0)
    @test bgc(1, 1, 1, grid, Val(:DOM), clock, zero_fields, auxiliary_fields) == 0.0
end
