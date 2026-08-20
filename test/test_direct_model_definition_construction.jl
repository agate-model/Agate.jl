using Oceananigans.Biogeochemistry:
    required_biogeochemical_auxiliary_fields, required_biogeochemical_tracers
using Test

using Agate.Configuration: AssimilationBinary, PalatabilityAllometric, Population, Pool
using Agate.Construction: construct
using Agate.Factories:
    DerivedDefault, FillDefault, ParameterDefinition, ParameterProvision, ParameterSpec
using Agate.Processes:
    Grazing, Growth, Light, ModelDefinition, NutrientResponse

function direct_npz_definition()
    provision(process, path, formulation, slot; qualifier=NamedTuple()) =
        ParameterProvision(process, path, formulation, slot; qualifier)
    definition(name, shape, provider, provides; axes=nothing) = ParameterDefinition(
        ParameterSpec(name, shape; axes, provides), provider
    )
    vector_definition(name, value, provides) =
        definition(name, :vector, FillDefault(value), provides; axes=:plankton)

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
            :preferential;
            consumer=:Z,
            resource=:P,
            unassimilated_destination=:N,
        ),
    )
    parameters = (
        vector_definition(
            :maximum_growth_rate,
            2e-5,
            provision(:growth_P, (:factors, :light), :smith, :maximum_rate),
        ),
        vector_definition(
            :alpha,
            2e-6,
            provision(:growth_P, (:factors, :light), :smith, :alpha),
        ),
        vector_definition(
            :nutrient_half_saturation,
            0.2,
            provision(
                :growth_P,
                (:factors, :nutrients),
                :monod,
                :K;
                qualifier=(resource=:N,),
            ),
        ),
        vector_definition(
            :maximum_predation_rate,
            5e-5,
            provision(:grazing_Z_on_P, (), :preferential, :maximum_rate),
        ),
        vector_definition(
            :holling_half_saturation,
            0.1,
            provision(:grazing_Z_on_P, (), :preferential, :half_saturation),
        ),
        vector_definition(
            :optimum_predator_prey_ratio,
            10.0,
            provision(
                :grazing_Z_on_P,
                (:palatability, :default),
                :allometric,
                :optimum_predator_prey_ratio,
            ),
        ),
        vector_definition(
            :specificity,
            0.3,
            provision(
                :grazing_Z_on_P, (:palatability, :default), :allometric, :specificity
            ),
        ),
        vector_definition(
            :protection,
            0.0,
            provision(
                :grazing_Z_on_P, (:palatability, :default), :allometric, :protection
            ),
        ),
        vector_definition(
            :assimilation_efficiency,
            0.7,
            provision(
                :grazing_Z_on_P,
                (:assimilation, :default),
                :binary,
                :assimilation_efficiency,
            ),
        ),
        definition(
            :palatability_matrix,
            :matrix,
            DerivedDefault(
                PalatabilityAllometric();
                deps=(:optimum_predator_prey_ratio, :specificity, :protection),
            ),
            provision(:grazing_Z_on_P, (), :preferential, :palatability);
            axes=(:consumer, :prey),
        ),
        definition(
            :assimilation_matrix,
            :matrix,
            DerivedDefault(AssimilationBinary(); deps=(:assimilation_efficiency,)),
            provision(:grazing_Z_on_P, (), :preferential, :assimilation);
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
