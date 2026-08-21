using Agate

const NiPiZD = Agate.Models.NiPiZD

using Agate.Library.Light: FunctionFieldPAR
using Agate.Diagnostics: box_model_mass_balance
using Agate.Configuration: Population, Pool
using Agate.Construction: construct
using Agate.Parameters:
    ConstDefault, FillDefault, ParameterDefinition, ParameterProvision
using Agate.Processes:
    FixedStoichiometry, Growth, Light, ModelDefinition, Nutrients, NutrientResponse

using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans.Units: day, minutes

function multi_nutrient_test_model(grid)
    components = (
        P=Population(; currency=:carbon, size_structure=[1.0]),
        DIC=Pool(:carbon),
        DIN=Pool(:nitrogen),
        PO4=Pool(:phosphorus),
    )
    processes = (
        growth_P=Growth(;
            population=:P,
            source=:DIC,
            factors=(
                light=Light(:geider; driver=:PAR),
                nutrients=Nutrients(
                    :liebig;
                    responses=(
                        nitrogen=NutrientResponse(:monod; resource=:DIN),
                        phosphorus=NutrientResponse(:monod; resource=:PO4),
                    ),
                ),
            ),
            stoichiometry=FixedStoichiometry(; reference=:carbon),
        ),
    )
    parameters = (
        ParameterDefinition(
            :nitrogen_to_carbon,
            ConstDefault(16 / 106);
            provides=ParameterProvision(:growth_P, :ratio; qualifier=(currency=:nitrogen,)),
        ),
        ParameterDefinition(
            :phosphorus_to_carbon,
            ConstDefault(1 / 106);
            provides=ParameterProvision(:growth_P, :ratio; qualifier=(currency=:phosphorus,)),
        ),
        ParameterDefinition(
            :maximum_growth_rate,
            FillDefault(2 / 86400);
            shape=:vector, axes=:plankton,
            provides=ParameterProvision(:growth_P, :maximum_rate),
        ),
        ParameterDefinition(
            :half_saturation_DIN,
            FillDefault(0.5);
            shape=:vector, axes=:plankton,
            provides=ParameterProvision(:growth_P, :K; qualifier=(resource=:DIN,)),
        ),
        ParameterDefinition(
            :half_saturation_PO4,
            FillDefault(0.5);
            shape=:vector, axes=:plankton,
            provides=ParameterProvision(:growth_P, :K; qualifier=(resource=:PO4,)),
        ),
        ParameterDefinition(
            :photosynthetic_slope,
            FillDefault(0.1 / 86400);
            shape=:vector, axes=:plankton,
            provides=ParameterProvision(:growth_P, :alpha),
        ),
        ParameterDefinition(
            :chlorophyll_to_carbon_ratio,
            FillDefault(0.02);
            shape=:vector, axes=:plankton,
            provides=ParameterProvision(:growth_P, :chlorophyll_to_carbon_ratio),
        ),
    )

    return construct(ModelDefinition(; components, processes, parameters); grid)
end

@testset "mass_balance" begin
    function build_box_model(bgc_instance, grid)
        bgc_model = Biogeochemistry(
            bgc_instance; light_attenuation=FunctionFieldPAR(; grid)
        )
        return BoxModel(; biogeochemistry=bgc_model, grid)
    end

    dt = 10minutes
    nsteps = Int(10day ÷ dt)

    @testset "NiPiZD box model" begin
        grid = BoxModelGrid()
        bgc_instance = NiPiZD.construct(; grid)
        box_model = build_box_model(bgc_instance, grid)
        set!(box_model; N=7, P_1=0.01, P_2=0.01, Z_1=0.05, Z_2=0.05, D=0.0)

        budgets = (total=[:N => 1, :P_1 => 1, :P_2 => 1, :Z_1 => 1, :Z_2 => 1, :D => 1],)

        result = box_model_mass_balance(box_model, budgets; dt, nsteps)
        @test isapprox(result.initial.total, result.final.total; rtol=1e-12, atol=0.0)
    end

    @testset "NiPiZD Float32 conservation" begin
        grid = BoxModelGrid(Float32)
        bgc_instance = NiPiZD.construct(; grid, scalar_type=Float32)
        box_model = build_box_model(bgc_instance, grid)
        set!(box_model; N=7f0, P_1=0.01f0, P_2=0.01f0, Z_1=0.05f0, Z_2=0.05f0, D=0f0)

        budgets = (total=[:N => 1, :P_1 => 1, :P_2 => 1, :Z_1 => 1, :Z_2 => 1, :D => 1],)

        result = box_model_mass_balance(box_model, budgets; dt, nsteps)
        @test result.initial.total isa Float32
        @test isapprox(result.initial.total, result.final.total; rtol=1e-5, atol=0.0)
    end

    @testset "Generic multi-nutrient conservation" begin
        grid = BoxModelGrid()
        bgc_instance = multi_nutrient_test_model(grid)
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
