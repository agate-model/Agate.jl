using Agate
using Agate.Library.Light: FunctionFieldPAR
using Agate.Introspection: tracer_names
using OceanBioME: Biogeochemistry, BoxModelGrid, BoxModel
using Oceananigans: set!, time_step!
using Oceananigans.Units: day
using Test

const NiPiZDReference = Agate.Models.NiPiZD
const HISTORICAL_NIPIZD_PARAMETERS = (
    mu=0.6989 / day,
    nutrient_half_saturation=2.3868,
    phytoplankton_mortality_to_nutrient=0.066 / day,
    zooplankton_mortality_to_nutrient=0.0102 / day,
    phytoplankton_mortality_to_detritus=0.0101 / day,
    maximum_grazing_rate=2.1522 / day,
    grazing_half_saturation=0.5573,
    assimilation_efficiency=0.9116,
    zooplankton_quadratic_mortality_to_detritus=0.3395 / day,
    detritus_remineralization=0.1213 / day,
    photosynthetic_slope=0.1953 / day,
)

reference_path = joinpath(@__DIR__, "reference", "nipizd_box_model.csv")
reference_rows = filter(
    row -> !isempty(row) && !startswith(row, '#') && row != "day,P",
    readlines(reference_path),
)
reference = map(reference_rows) do row
    fields = split(row, ',')
    return (day=parse(Int, fields[1]), P=parse(Float64, fields[2]))
end

function historical_nipizd_bgc()
    p = HISTORICAL_NIPIZD_PARAMETERS
    return NiPiZDReference.construct(;
        size_structure=(
            phytoplankton=(P=[1.0],),
            zooplankton=(Z=[1.0],),
        ),
        parameters=(
            detritus_remineralization=p.detritus_remineralization,
            linear_mortality=(
                P_1=p.phytoplankton_mortality_to_nutrient,
                Z_1=p.zooplankton_mortality_to_nutrient,
            ),
            linear_detrital_mortality=(P_1=p.phytoplankton_mortality_to_detritus,),
            quadratic_mortality=(Z_1=p.zooplankton_quadratic_mortality_to_detritus,),
            maximum_growth_rate=(P_1=p.mu,),
            nutrient_half_saturation=(P_1=p.nutrient_half_saturation,),
            alpha=(P_1=p.photosynthetic_slope,),
            maximum_predation_rate=(Z_1=p.maximum_grazing_rate,),
            holling_half_saturation=(Z_1=p.grazing_half_saturation,),
        ),
        palatability_matrix=[1.0;;],
        assimilation_matrix=[p.assimilation_efficiency;;],
    )
end

@testset "NiPiZD historical scientific trajectory" begin
    @test [entry.day for entry in reference] == collect(10:10:1000)

    grid = BoxModelGrid()
    bgc = historical_nipizd_bgc()
    @test tracer_names(bgc) == [:N, :D, :Z_1, :P_1]

    bgc_model = Biogeochemistry(
        bgc; light_attenuation=FunctionFieldPAR(; grid)
    )
    box_model = BoxModel(; biogeochemistry=bgc_model)
    set!(box_model; N=7.0, P_1=0.01, Z_1=0.05, D=0.0)

    for i in 1:1000
        time_step!(box_model, 1day)
        i % 10 == 0 || continue
        P = box_model.fields.P_1.data[1, 1, 1]
        @test P ≈ reference[i ÷ 10].P rtol=1e-12 atol=0
    end

    total_N = sum(name -> box_model.fields[name].data[1, 1, 1], (:N, :P_1, :Z_1, :D))
    @test total_N ≈ 7.06 rtol=1e-12 atol=1e-12
end
