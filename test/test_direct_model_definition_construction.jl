using Oceananigans.Biogeochemistry:
    required_biogeochemical_auxiliary_fields, required_biogeochemical_tracers
using Test

using Agate.Components: Plankton, Pool
using Agate.Parameters: AllometricPalatability, ConsumerAssimilation, ConstantDefault, ConstructionParameter, DerivedDefault, Parameter
using Agate.Construction: construct
using Agate.Introspection: plankton_diameters
using Agate.Processes:
    Consumption, FixedStoichiometry, Growth, Light, LinearMortality, ModelDefinition, Mortality,
    NutrientResponse, Products, Smith, Monod, PreferentialGrazing

function direct_npz_definition()
    components = (
        N=Pool(:nitrogen),
        P=Plankton(; states=(nitrogen=:nitrogen,), reference_state=:nitrogen, size_structure=[1.0]),
        Z=Plankton(; states=(nitrogen=:nitrogen,), reference_state=:nitrogen, size_structure=[10.0]),
    )
    processes = (
        growth_P=Growth(;
            plankton=:P,
            reference_resource=:N,
            bindings=(maximum_rate=:maximum_growth_rate,),
            factors=(
                light=Light(Smith(); driver=:PAR),
                nutrients=NutrientResponse(
                    Monod(); resource=:N,
                    bindings=(half_saturation=:nutrient_half_saturation,)
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
        optimum_predator_prey_ratio=ConstructionParameter(
            ConstantDefault(10.0);
            axes=:plankton,
        ),
        specificity=ConstructionParameter(
            ConstantDefault(0.3);
            axes=:plankton,
        ),
        protection=ConstructionParameter(
            ConstantDefault(0.0);
            axes=:plankton,
        ),
        assimilation_efficiency=ConstructionParameter(
            ConstantDefault(0.7);
            axes=:plankton,
        ),
        palatability_matrix=Parameter(
            DerivedDefault(
                AllometricPalatability();
                deps=(:optimum_predator_prey_ratio, :specificity, :protection),
            )
        ),
        assimilation_matrix=Parameter(
            DerivedDefault(ConsumerAssimilation(); deps=(:assimilation_efficiency,))
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
    tendencies = values(model_tendencies(bgc, args))
    @test all(isfinite, tendencies)
    @test isapprox(sum(tendencies), 0; atol=10 * eps(sum(abs, tendencies)))
end

@testset "Process-local interaction defaults and modeled edges" begin
    components = (
        P1=Plankton(; states=(nitrogen=:nitrogen,), reference_state=:nitrogen, size_structure=[1.0]),
        P2=Plankton(; states=(nitrogen=:nitrogen,), reference_state=:nitrogen, size_structure=[4.0]),
        Z1=Plankton(; states=(nitrogen=:nitrogen,), reference_state=:nitrogen, size_structure=[10.0]),
        Z2=Plankton(; states=(nitrogen=:nitrogen,), reference_state=:nitrogen, size_structure=[20.0]),
    )
    common = (maximum_rate=:maximum_predation_rate, half_saturation=:half_saturation, palatability=:palatability_shared)
    processes = (
        a=Consumption(PreferentialGrazing(); consumers=:Z1, resources=:P1, bindings=(; common..., assimilation=:assimilation_other)),
        b=Consumption(PreferentialGrazing(); consumers=:Z2, resources=:P2, bindings=(; common..., assimilation=:assimilation_local)),
    )
    parameters = (
        maximum_predation_rate=Parameter(1.0), half_saturation=Parameter(1.0),
        palatability_shared=Parameter(0.5), assimilation_other=Parameter(0.5),
        assimilation_local=Parameter(DerivedDefault(ConsumerAssimilation(); deps=(:assimilation_efficiency,))),
        assimilation_efficiency=ConstructionParameter(0.5; axes=:plankton),
    )
    bgc = construct(
        ModelDefinition(; components, processes, parameters);
        parameter_overrides=(assimilation_efficiency=[0.0, 0.0, 0.6, 0.7],),
    )

    @test bgc.parameters.assimilation_local == fill(0.7, 1, 1)
    shared = Agate.Introspection.interaction_matrix(bgc, :palatability_shared)
    @test Set(shared.edges) == Set(((:Z1_1, :P1_1), (:Z2_1, :P2_1)))
    @test_throws ArgumentError Agate.Runtime.active_parameters(
        bgc; palatability_shared=((:Z1_1, :P2_1),)
    )

end

@testset "Multi-element products" begin
    components = (
        P=Plankton(; states=(carbon=:carbon,), reference_state=:carbon, size_structure=[1.0]),
        DOC=Pool(:carbon), POC=Pool(:carbon), XOC=Pool(:carbon),
        DON=Pool(:nitrogen), PON=Pool(:nitrogen), XON=Pool(:nitrogen),
    )
    stoichiometry = FixedStoichiometry(;
        reference_element=:carbon,
        bindings=(ratio=(nitrogen=:nitrogen_to_carbon,),),
    )
    processes = (
        mortality_P=Mortality(
            LinearMortality();
            plankton=:P,
            bindings=(rate=:linear_mortality,),
            products=Products(
                (
                    DOM=(carbon=:DOC, nitrogen=:DON),
                    POM=(carbon=:POC, nitrogen=:PON),
                    XOM=(carbon=:XOC, nitrogen=:XON),
                );
                fractions=(POM=:POM_fraction, XOM=:XOM_fraction),
                stoichiometry,
            ),
        ),
    )
    parameters = (
        linear_mortality=Parameter(ConstantDefault(1e-6)),
        POM_fraction=Parameter(ConstantDefault(0.25)),
        XOM_fraction=Parameter(ConstantDefault(0.25)),
        nitrogen_to_carbon=Parameter(ConstantDefault(0.2)),
    )
    bgc = construct(ModelDefinition(; components, processes, parameters); grid=dummy_grid(Float64))

    tracers = required_biogeochemical_tracers(bgc)
    state = Dict(tracer => (tracer === :P_1 ? 1.0 : 0.0) for tracer in tracers)
    args = (0.0, 0.0, 0.0, 0.0, (state[tracer] for tracer in tracers)...)
    tendency = model_tendencies(bgc, args; tracers)
    actual = (tendency.DOC, tendency.POC, tendency.XOC,
              tendency.DON, tendency.PON, tendency.XON, tendency.P_1)
    expected = (0.5e-6, 0.25e-6, 0.25e-6,
                0.1e-6, 0.05e-6, 0.05e-6, -1e-6)
    @test all(process_compiler_isapprox.(actual, expected))
    @test_throws ArgumentError construct(
        ModelDefinition(; components, processes, parameters); parameter_overrides=(POM_fraction=0.8, XOM_fraction=0.8))

    active = Agate.Runtime.active_parameters(bgc; POM_fraction=true, XOM_fraction=true)
    @test_throws ArgumentError Agate.Runtime.parameterized(
        bgc, [0.8, 0.8]; active_parameters=active)
end
