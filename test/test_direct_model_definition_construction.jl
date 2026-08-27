using Oceananigans.Biogeochemistry:
    required_biogeochemical_auxiliary_fields, required_biogeochemical_tracers
using Test

using Agate.Configuration: AssimilationBinary, PalatabilityAllometric, Population, Pool
using Agate.Construction: construct
using Agate.Introspection: plankton_diameters
using Agate.Parameters: DerivedDefault, ConstantDefault, MetaParameter, Parameter
using Agate.Processes:
    Consumption, FixedStoichiometry, Growth, Light, LinearMortality, ModelDefinition, Mortality,
    NutrientResponse, Products, Smith, Monod, PreferentialGrazing

function direct_npz_definition()
    components = (
        N=Pool(:nitrogen),
        P=Population(:nitrogen; size_structure=[1.0]),
        Z=Population(:nitrogen; size_structure=[10.0]),
    )
    processes = (
        growth_P=Growth(;
            populations=:P,
            bindings=(maximum_rate=:maximum_growth_rate,),
            factors=(
                light=Light(Smith(); driver=:PAR),
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
        maximum_growth_rate=Parameter(ConstantDefault(2e-5)),
        alpha=Parameter(ConstantDefault(2e-6)),
        nutrient_half_saturation=Parameter(ConstantDefault(0.2)),
        maximum_predation_rate=Parameter(ConstantDefault(5e-5)),
        holling_half_saturation=Parameter(ConstantDefault(0.1)),
        optimum_predator_prey_ratio=MetaParameter(
            ConstantDefault(10.0);
            axes=:plankton,
        ),
        specificity=MetaParameter(
            ConstantDefault(0.3);
            axes=:plankton,
        ),
        protection=MetaParameter(
            ConstantDefault(0.0);
            axes=:plankton,
        ),
        assimilation_efficiency=MetaParameter(
            ConstantDefault(0.7);
            axes=:plankton,
        ),
        palatability_matrix=Parameter(
            DerivedDefault(
                PalatabilityAllometric();
                deps=(:optimum_predator_prey_ratio, :specificity, :protection),
            )
        ),
        assimilation_matrix=Parameter(
            DerivedDefault(AssimilationBinary(); deps=(:assimilation_efficiency,))
        ),
    )
    return ModelDefinition(; components, processes, parameters)
end

@testset "Direct construction from ModelDefinition" begin
    definition = direct_npz_definition()
    bgc = construct(
        definition;
        grid=dummy_grid(Float64),
        parameter_overrides=(maximum_growth_rate=[3e-5],),
    )

    @test required_biogeochemical_tracers(bgc) == (:N, :P_1, :Z_1)
    @test required_biogeochemical_auxiliary_fields(bgc) == (:PAR,)
    @test bgc.parameters.maximum_growth_rate == [3e-5]
    @test bgc.parameters.alpha == [2e-6]
    @test bgc.parameters.maximum_predation_rate == [5e-5]
    @test plankton_diameters(bgc) == [1.0, 10.0]
    @test size(bgc.parameters.palatability_matrix) ==
          size(bgc.parameters.assimilation_matrix) == (1, 1)
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
        linear_mortality=Parameter(ConstantDefault(1e-6)),
        fraction_a=Parameter(ConstantDefault(0.4)),
        fraction_b=Parameter(ConstantDefault(nextfloat(0.6))),
    )
    definition = ModelDefinition(; components, processes, parameters)

    bgc = construct(definition; grid=dummy_grid(Float64))
    @test hasproperty(bgc.parameters, :fraction_a)
    @test !hasproperty(bgc.parameters, :fraction_b)
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

@testset "Multi-currency products" begin
    components = (
        P=Population(:carbon; size_structure=[1.0]),
        DOC=Pool(:carbon),
        POC=Pool(:carbon),
        DON=Pool(:nitrogen),
        PON=Pool(:nitrogen),
    )
    stoichiometry = FixedStoichiometry(;
        reference=:carbon,
        bindings=(ratio=(nitrogen=:nitrogen_to_carbon,),),
    )
    processes = (
        mortality_P=Mortality(
            LinearMortality();
            populations=:P,
            bindings=(rate=:linear_mortality,),
            products=Products(
                (
                    DOM=(carbon=:DOC, nitrogen=:DON),
                    POM=(carbon=:POC, nitrogen=:PON),
                );
                fractions=(POM=:POM_fraction,),
                stoichiometry,
            ),
        ),
    )
    parameters = (
        linear_mortality=Parameter(ConstantDefault(1e-6)),
        POM_fraction=Parameter(ConstantDefault(0.25)),
        nitrogen_to_carbon=Parameter(ConstantDefault(0.2)),
    )
    bgc = construct(ModelDefinition(; components, processes, parameters); grid=dummy_grid(Float64))

    tracers = required_biogeochemical_tracers(bgc)
    state = Dict(:DOC => 0.0, :POC => 0.0, :DON => 0.0, :PON => 0.0, :P_1 => 1.0)
    args = (0.0, 0.0, 0.0, 0.0, (state[tracer] for tracer in tracers)...)
    tendency = NamedTuple{tracers}(Tuple(bgc(Val(tracer), args...) for tracer in tracers))

    actual = (tendency.DOC, tendency.POC, tendency.DON, tendency.PON, tendency.P_1)
    expected = (0.75e-6, 0.25e-6, (0.75 * 0.2) * 1e-6, (0.25 * 0.2) * 1e-6, -1e-6)
    @test all(process_compiler_isapprox.(actual, expected))
end
