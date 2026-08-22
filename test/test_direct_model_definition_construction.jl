using Oceananigans.Biogeochemistry:
    required_biogeochemical_auxiliary_fields, required_biogeochemical_tracers
using Test

using Agate.Configuration: AssimilationBinary, PalatabilityAllometric, Population, Pool
using Agate.Construction: construct
using Agate.Parameters:
    DerivedDefault, ConstantDefault, Parameter, ParameterProvision
using Agate.Processes:
    Consumption, Growth, Light, ModelDefinition, NutrientResponse, ProductRouting,
    Smith, Monod, PreferentialGrazing, DirectRouting

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
                light=Light(Smith(); driver=:PAR),
                nutrients=NutrientResponse(Monod(); resource=:N),
            ),
        ),
        grazing_Z_on_P=Consumption(
            PreferentialGrazing();
            consumers=:Z,
            resources=:P,
            routing=ProductRouting(DirectRouting(); destination=:N),
        ),
    )
    parameters = (
        maximum_growth_rate=Parameter(
            ConstantDefault(2e-5);
            axes=:plankton,
            provides=ParameterProvision(:growth_P, :maximum_rate),
        ),
        alpha=Parameter(
            ConstantDefault(2e-6);
            axes=:plankton,
            provides=ParameterProvision(:growth_P, :alpha),
        ),
        nutrient_half_saturation=Parameter(
            ConstantDefault(0.2);
            axes=:plankton,
            provides=ParameterProvision(:growth_P, :K),
        ),
        maximum_predation_rate=Parameter(
            ConstantDefault(5e-5);
            axes=:plankton,
            provides=ParameterProvision(:grazing_Z_on_P, :maximum_rate),
        ),
        holling_half_saturation=Parameter(
            ConstantDefault(0.1);
            axes=:plankton,
            provides=ParameterProvision(:grazing_Z_on_P, :half_saturation),
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
            provides=ParameterProvision(:grazing_Z_on_P, :palatability),
        ),
        assimilation_matrix=Parameter(
            DerivedDefault(AssimilationBinary(); deps=(:assimilation_efficiency,));
            axes=(:consumer, :prey),
            provides=ParameterProvision(:grazing_Z_on_P, :assimilation),
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
