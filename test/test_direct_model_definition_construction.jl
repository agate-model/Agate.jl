using Oceananigans.Biogeochemistry:
    required_biogeochemical_auxiliary_fields, required_biogeochemical_tracers
using Test

using Agate.Configuration: AssimilationBinary, PalatabilityAllometric, Population, Pool
using Agate.Construction: construct
using Agate.Parameters:
    DerivedDefault, ConstantDefault, ParameterDefinition, ParameterProvision
using Agate.Processes:
    Grazing, Growth, Light, ModelDefinition, NutrientResponse

function direct_npz_definition()
    components = (
        N=Pool(:nitrogen),
        P=Population(; currency=:nitrogen, size_structure=[1.0]),
        Z=Population(; currency=:nitrogen, size_structure=[10.0]),
    )
    processes = (
        growth_P=Growth(;
            population=:P,
            factors=(
                light=Light(:smith; driver=:PAR),
                nutrients=NutrientResponse(:monod; resource=:N),
            ),
        ),
        grazing_Z_on_P=Grazing(
            :preferential; consumer=:Z, resource=:P, unassimilated_destination=:N
        ),
    )
    parameters = (
        ParameterDefinition(
            :maximum_growth_rate,
            ConstantDefault(2e-5);
            axes=:plankton,
            provides=ParameterProvision(:growth_P, :maximum_rate),
        ),
        ParameterDefinition(
            :alpha,
            ConstantDefault(2e-6);
            axes=:plankton,
            provides=ParameterProvision(:growth_P, :alpha),
        ),
        ParameterDefinition(
            :nutrient_half_saturation,
            ConstantDefault(0.2);
            axes=:plankton,
            provides=ParameterProvision(:growth_P, :K),
        ),
        ParameterDefinition(
            :maximum_predation_rate,
            ConstantDefault(5e-5);
            axes=:plankton,
            provides=ParameterProvision(:grazing_Z_on_P, :maximum_rate),
        ),
        ParameterDefinition(
            :holling_half_saturation,
            ConstantDefault(0.1);
            axes=:plankton,
            provides=ParameterProvision(:grazing_Z_on_P, :half_saturation),
        ),
        ParameterDefinition(
            :optimum_predator_prey_ratio,
            ConstantDefault(10.0);
            axes=:plankton,
        ),
        ParameterDefinition(
            :specificity,
            ConstantDefault(0.3);
            axes=:plankton,
        ),
        ParameterDefinition(
            :protection,
            ConstantDefault(0.0);
            axes=:plankton,
        ),
        ParameterDefinition(
            :assimilation_efficiency,
            ConstantDefault(0.7);
            axes=:plankton,
        ),
        ParameterDefinition(
            :palatability_matrix,
            DerivedDefault(
                PalatabilityAllometric();
                deps=(:optimum_predator_prey_ratio, :specificity, :protection),
            );
            axes=(:consumer, :prey),
            runtime_path=(:interactions, :palatability),
            provides=ParameterProvision(:grazing_Z_on_P, :palatability),
        ),
        ParameterDefinition(
            :assimilation_matrix,
            DerivedDefault(AssimilationBinary(); deps=(:assimilation_efficiency,));
            axes=(:consumer, :prey),
            runtime_path=(:interactions, :assimilation),
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
