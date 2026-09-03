using Agate

using Agate.Diagnostics: box_model_mass_balance
using Agate.Components: Plankton, Pool
using Agate.Construction: construct
using Agate.Parameters: ConstantDefault, Parameter
using Agate.Processes:
    FixedStoichiometry, Growth, Light, ModelDefinition, NutrientLimitation, NutrientResponse,
    Geider, Liebig, FrankTNorm, Monod

using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans.Fields: ConstantField
using Oceananigans.Units: day, minutes

function multi_nutrient_test_model(grid; nutrient_formulation=Liebig())
    components = (
        P=Plankton(; states=(carbon=:carbon,), reference_state=:carbon, size_structure=[1.0]),
        DIC=Pool(:carbon),
        DIN=Pool(:nitrogen),
        PO4=Pool(:phosphorus),
    )
    processes = (
        growth_P=Growth(;
            plankton=:P,
            reference_resource=:DIC,
            additional_resources=(nitrogen=:DIN, phosphorus=:PO4),
            bindings=(maximum_rate=:maximum_growth_rate,),
            factors=(
                light=Light(
                    Geider();
                    driver=:PAR,
                    bindings=(alpha=:photosynthetic_slope,),
                ),
                nutrients=NutrientLimitation(
                    nutrient_formulation;
                    bindings=nutrient_formulation isa FrankTNorm ?
                             (sharpness=:frank_sharpness,) : NamedTuple(),
                    responses=(
                        nitrogen=NutrientResponse(
                            Monod(); resource=:DIN,
                            bindings=(half_saturation=:half_saturation_DIN,)
                        ),
                        phosphorus=NutrientResponse(
                            Monod(); resource=:PO4,
                            bindings=(half_saturation=:half_saturation_PO4,)
                        ),
                    ),
                ),
            ),
            stoichiometry=FixedStoichiometry(;
                reference_element=:carbon,
                bindings=(
                    ratio=(
                        nitrogen=:nitrogen_to_carbon,
                        phosphorus=:phosphorus_to_carbon,
                    ),
                ),
            ),
        ),
    )
    frank_tnorm_parameters = nutrient_formulation isa FrankTNorm ? (
        frank_sharpness=Parameter(ConstantDefault(25)),
    ) : (;)
    parameters = merge((
        nitrogen_to_carbon=Parameter(ConstantDefault(16 / 106)),
        phosphorus_to_carbon=Parameter(ConstantDefault(1 / 106)),
        maximum_growth_rate=Parameter(ConstantDefault(2 / 86400)),
        half_saturation_DIN=Parameter(ConstantDefault(0.5)),
        half_saturation_PO4=Parameter(ConstantDefault(0.5)),
        photosynthetic_slope=Parameter(ConstantDefault(0.1 / 86400)),
        chlorophyll_to_carbon_ratio=Parameter(ConstantDefault(0.02)),
    ), frank_tnorm_parameters)

    return construct(ModelDefinition(; components, processes, parameters); grid)
end

@testset "mass_balance" begin
    function build_box_model(bgc_instance, grid)
        light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(ConstantField(100.0))
        bgc_model = Biogeochemistry(bgc_instance; light_attenuation)
        return BoxModel(; biogeochemistry=bgc_model, grid)
    end

    dt = 10minutes
    nsteps = Int(10day ÷ dt)

    budgets = (total=[:N => 1, :P_1 => 1, :P_2 => 1, :Z_1 => 1, :Z_2 => 1, :D => 1],)
    for (T, rtol) in ((Float64, 1e-12), (Float32, 1e-5))
        @testset "NiPiZD $(T) conservation" begin
            grid = BoxModelGrid(T)
            bgc = Agate.Models.NiPiZD.construct(; grid, scalar_type=T)
            box_model = build_box_model(bgc, grid)
            set!(box_model; N=T(7), P_1=T(0.01), P_2=T(0.01),
                 Z_1=T(0.05), Z_2=T(0.05), D=zero(T))

            result = box_model_mass_balance(box_model, budgets; dt, nsteps)
            @test result.initial.total isa T
            @test isapprox(result.initial.total, result.final.total; rtol, atol=0.0)
        end
    end

    @testset "Quota model elemental conservation" begin
        grid = BoxModelGrid()
        bgc_instance = construct(quota_definition(); grid)
        box_model = build_box_model(bgc_instance, grid)
        set!(
            box_model;
            DIC=10.0,
            DIN=1.0,
            PO4=1.0,
            P_1_carbon=1.0,
            P_1_nitrogen=0.1,
            P_1_phosphorus=0.01,
            P_2_carbon=0.5,
            P_2_nitrogen=0.075,
            P_2_phosphorus=0.0075,
        )

        quota_budgets = (
            carbon=[:DIC => 1, :P_1_carbon => 1, :P_2_carbon => 1],
            nitrogen=[:DIN => 1, :P_1_nitrogen => 1, :P_2_nitrogen => 1],
            phosphorus=[:PO4 => 1, :P_1_phosphorus => 1, :P_2_phosphorus => 1],
        )

        # The shared quota fixture uses rates in inverse seconds, so integrate over a
        # short interval rather than the multi-day NiPiZD conservation horizon.
        result = box_model_mass_balance(box_model, quota_budgets; dt=0.01, nsteps=100)
        @test isapprox(result.initial.carbon, result.final.carbon; rtol=1e-10)
        @test isapprox(result.initial.nitrogen, result.final.nitrogen; rtol=1e-10)
        @test isapprox(result.initial.phosphorus, result.final.phosphorus; rtol=1e-10)
    end

    @testset "Generic multi-nutrient conservation" begin
        for nutrient_formulation in (Liebig(), FrankTNorm())
            @testset "$(nameof(typeof(nutrient_formulation)))" begin
                grid = BoxModelGrid()
                bgc_instance = multi_nutrient_test_model(grid; nutrient_formulation)
                if nutrient_formulation isa FrankTNorm
                    @test bgc_instance.parameters.frank_sharpness == 25
                    state = (DIC=2.0, DIN=7.0, PO4=3.0, P_1=0.01)
                    tracers = Tuple(
                        getproperty(state, name)
                        for name in Agate.Introspection.tracer_names(bgc_instance)
                    )
                    args = (0, 0, 0, 0, tracers..., 100.0)
                    limitations = (
                        Agate.Processes.factor_value(Monod(), state.DIN, 0.5),
                        Agate.Processes.factor_value(Monod(), state.PO4, 0.5),
                    )
                    light = Agate.Processes.factor_value(
                        Geider(), 100.0, 2 / 86400, 0.1 / 86400, 0.02
                    )
                    nutrients = Agate.Processes.factor_value(FrankTNorm(), limitations, 25)
                    @test bgc_instance(Val(:P_1), args...) ≈
                          state.P_1 * (2 / 86400) * light * nutrients
                end
                box_model = build_box_model(bgc_instance, grid)
                set!(box_model; DIC=2.0, DIN=7.0, PO4=3.0, P_1=0.01)

                n2c = bgc_instance.parameters.nitrogen_to_carbon
                p2c = bgc_instance.parameters.phosphorus_to_carbon
                budgets = (
                    carbon=[:DIC => 1, :P_1 => 1],
                    nitrogen=[:DIN => 1, :P_1 => n2c],
                    phosphorus=[:PO4 => 1, :P_1 => p2c],
                )

                result = box_model_mass_balance(box_model, budgets; dt, nsteps)
                @test isapprox(result.initial.carbon, result.final.carbon; rtol=1e-6)
                @test isapprox(result.initial.nitrogen, result.final.nitrogen; rtol=1e-6)
                @test isapprox(result.initial.phosphorus, result.final.phosphorus; rtol=1e-6)
            end
        end
    end

end
