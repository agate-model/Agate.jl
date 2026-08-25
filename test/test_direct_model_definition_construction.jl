using Oceananigans.Biogeochemistry:
    required_biogeochemical_auxiliary_fields, required_biogeochemical_tracers
using Test

using Agate.Configuration: AssimilationBinary, PalatabilityAllometric, Population, Pool
using Agate.Construction: construct
using Agate.Parameters: DerivedDefault, ConstantDefault, Parameter
using Agate.Processes:
    Consumption, Growth, Light, LinearMortality, ModelDefinition, Mortality, NutrientResponse,
    Products, Smith, Monod, PreferentialGrazing

function direct_npz_definition()
    components = (
        N=Pool(:nitrogen),
        P=Population(:nitrogen; size_structure=[1.0]),
        Z=Population(:nitrogen; size_structure=[10.0]),
    )
    processes = (
        growth_P=Growth(;
            populations=:P,
            factors=(
                light=Light(
                    Smith(); driver=:PAR, bindings=(maximum_rate=:maximum_growth_rate,)
                ),
                nutrients=NutrientResponse(
                    Monod(); resource=:N, bindings=(K=:nutrient_half_saturation,)
                ),
            ),
        ),
        grazing_Z_on_P=Consumption(
            PreferentialGrazing();
            consumers=:Z,
            resources=:P,
            bindings=(
                maximum_rate=:maximum_predation_rate,
                half_saturation=:holling_half_saturation,
                palatability=:palatability_matrix,
                assimilation=:assimilation_matrix,
            ),
            unassimilated_products=:N,
        ),
    )
    parameters = (
        maximum_growth_rate=Parameter(
            ConstantDefault(2e-5);
            axes=:plankton,
        ),
        alpha=Parameter(
            ConstantDefault(2e-6);
            axes=:plankton,
        ),
        nutrient_half_saturation=Parameter(
            ConstantDefault(0.2);
            axes=:plankton,
        ),
        maximum_predation_rate=Parameter(
            ConstantDefault(5e-5);
            axes=:plankton,
        ),
        holling_half_saturation=Parameter(
            ConstantDefault(0.1);
            axes=:plankton,
        ),
        optimum_predator_prey_ratio=Parameter(
            ConstantDefault(10.0);
            axes=:plankton,
        ),
        specificity=Parameter(
            ConstantDefault(0.3);
            axes=:plankton,
        ),
        protection=Parameter(
            ConstantDefault(0.0);
            axes=:plankton,
        ),
        assimilation_efficiency=Parameter(
            ConstantDefault(0.7);
            axes=:plankton,
        ),
        palatability_matrix=Parameter(
            DerivedDefault(
                PalatabilityAllometric();
                deps=(:optimum_predator_prey_ratio, :specificity, :protection),
            );
            axes=(:consumer, :prey),
        ),
        assimilation_matrix=Parameter(
            DerivedDefault(AssimilationBinary(); deps=(:assimilation_efficiency,));
            axes=(:consumer, :prey),
        ),
    )
    return ModelDefinition(; components, processes, parameters)
end

@testset "Direct construction from ModelDefinition" begin
    definition = direct_npz_definition()
    bgc = construct(
        definition;
        grid=dummy_grid(Float64),
        parameter_overrides=(maximum_growth_rate=[3e-5, 0.0],),
    )

    @test required_biogeochemical_tracers(bgc) == (:N, :P_1, :Z_1)
    @test required_biogeochemical_auxiliary_fields(bgc) == (:PAR,)
    @test bgc.parameters.maximum_growth_rate == [3e-5, 0.0]
    @test bgc.parameters.alpha == [2e-6, 0.0]
    @test bgc.parameters.maximum_predation_rate == [0.0, 5e-5]
    @test bgc.plankton_diameters == (1.0, 10.0)
    @test size(bgc.parameters.palatability_matrix) ==
          size(bgc.parameters.assimilation_matrix) == (1, 1)
    @test all(
        equation -> isbitstype(typeof(equation)) &&
                    all(term -> isbitstype(typeof(term)), equation.terms),
        values(bgc.tracer_functions),
    )

    args = (0.0, 0.0, 0.0, 0.0, 1.0, 0.1, 0.05, 100.0)
    tendencies = map(
        tracer -> bgc(Val(tracer), args...), required_biogeochemical_tracers(bgc)
    )
    @test all(isfinite, tendencies)
    @test isapprox(sum(tendencies), 0; atol=10 * eps(sum(abs, tendencies)))
end

@testset "Explicit product fractions" begin
    components = (
        P=Population(:nitrogen; size_structure=[1.0]),
        A=Pool(:nitrogen),
        B=Pool(:nitrogen),
    )
    processes = (
        mortality_P=Mortality(
            LinearMortality();
            populations=:P,
            bindings=(rate=:linear_mortality,),
            products=Products(
                (a=:A, b=:B);
                fractions=(a=:fraction_a, b=:fraction_b),
            ),
        ),
    )
    parameters = (
        linear_mortality=Parameter(ConstantDefault(1e-6); axes=:plankton),
        fraction_a=Parameter(ConstantDefault(0.4)),
        fraction_b=Parameter(ConstantDefault(nextfloat(0.6))),
    )
    definition = ModelDefinition(; components, processes, parameters)

    bgc = construct(definition; grid=dummy_grid(Float64))
    args = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0)
    tendencies = map(tracer -> bgc(Val(tracer), args...), (:A, :B, :P_1))
    conservative_b = (1.0 - 0.4) * 1e-6
    @test tendencies[2] == conservative_b
    @test tendencies[2] != nextfloat(0.6) * 1e-6
    @test isapprox(sum(tendencies), 0; atol=10 * eps(sum(abs, tendencies)))
    @test_throws ArgumentError construct(
        definition;
        grid=dummy_grid(Float64),
        parameter_overrides=(fraction_b=0.5,),
    )
end
