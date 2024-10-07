using Agate
using OceanBioME
using Oceananigans
using Oceananigans: Clock
using Oceananigans.Units
using Oceananigans.Fields: FunctionField

const year = years = 365day
const z = -10 # specify the nominal depth of the box for the PAR profile

@testset "NPZD box model" begin

    # ==================================================
    # Agate NPZD model
    # ==================================================

    # get parameter and tracer definitions from examples
    # the values here are the same as OceanBioME defaults
    include("../examples/NPZD/model_definition.jl")

    NPZD = create_bgc_struct(:NPZD, parameters)
    add_bgc_methods(
        NPZD,
        tracers,
        auxiliary_fields=aux_field_vars,
        helper_functions="../examples/NPZD/functions.jl"
    )
    npzd_model = NPZD()
    init_conditions = (N = 7.0, P = 0.01, Z = 0.05, D=0.0)
    agate_box_model = create_box_model(npzd_model, init_conditions)

    # ==================================================
    # OceanBioME NPZD model
    # ==================================================

    PAR⁰(t) = 60 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2
    PAR_f(t) = PAR⁰(t) * exp(0.2z) # Modify the PAR based on the nominal depth and exponential decay

    grid = BoxModelGrid()
    clock = Clock(time = zero(grid))
    PAR = FunctionField{Center, Center, Center}(PAR_f, grid; clock)

    biogeochemistry = NutrientPhytoplanktonZooplanktonDetritus(;
        grid,
        light_attenuation_model = PrescribedPhotosyntheticallyActiveRadiation(PAR),
        # this is probably not necessary but ensuring consistency here
        sinking_speeds = NamedTuple()
    )
    oceanbiome_box_model = BoxModel(;biogeochemistry, clock)
    set!(oceanbiome_box_model, N=7, P=0.01, Z=0.05, D=0.0)

    # ==================================================
    # Compare
    # ==================================================

    Δt = 1day
    for i in range(1,1000)
        time_step!(oceanbiome_box_model, Δt)
        time_step!(agate_box_model, Δt)
        if mod(i, 10) == 0
            @test agate_box_model.fields.P.data[1,1,1] === oceanbiome_box_model.fields.P.data[1,1,1]
        end
    end
end
