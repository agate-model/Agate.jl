using Agate
using Agate.Light
using OceanBioME
using Oceananigans
using Oceananigans: Clock
using Oceananigans.Units
using Oceananigans.Fields: FunctionField

const year = years = 365day

@testset "NPZD box model" begin

    # ==================================================
    # Agate NPZD model
    # ==================================================

    include("../examples/NPZD/model.jl")
    npzd_model = NPZD()
    init_conditions = (N = 7.0, P = 0.01, Z = 0.05, D = 0.0)
    agate_box_model = create_box_model(npzd_model, init_conditions)

    # ==================================================
    # OceanBioME NPZD model
    # ==================================================
    grid = BoxModelGrid()
    clock = Clock(time = zero(grid))
    PAR = FunctionField{Center,Center,Center}(PAR_box, grid; clock)

    biogeochemistry = NutrientPhytoplanktonZooplanktonDetritus(;
        grid,
        light_attenuation_model = PrescribedPhotosyntheticallyActiveRadiation(PAR),
        # this is probably not necessary but ensuring consistency here
        sinking_speeds = NamedTuple(),
    )
    oceanbiome_box_model = BoxModel(; biogeochemistry, clock)
    set!(oceanbiome_box_model, N = 7, P = 0.01, Z = 0.05, D = 0.0)

    # ==================================================
    # Compare
    # ==================================================

    Δt = 1day
    for i in range(1, 1000)
        time_step!(oceanbiome_box_model, Δt)
        time_step!(agate_box_model, Δt)
        if mod(i, 10) == 0
            @test agate_box_model.fields.P.data[1, 1, 1] ===
                  oceanbiome_box_model.fields.P.data[1, 1, 1]
        end
    end
end