using Test

using Agate.Configuration: Plankton, Pool
using Agate.Parameters: Parameter
using Agate.Processes:
    Consumption, Growth, LinearMortality, Monod, Mortality, ModelDefinition,
    NutrientResponse, PreferentialGrazing, Products

function _test_tendencies(model, state::NamedTuple, expected::NamedTuple)
    names = Agate.Introspection.tracer_names(model)
    args = (0.0, 0.0, 0.0, 0.0, Tuple(getproperty(state, name) for name in names)...)
    for (tracer, tendency) in pairs(expected)
        @test model(Val(tracer), args...) ≈ tendency
    end
end

@testset "Multi-state process semantics" begin
    @testset "non-quota growth permits physiological states" begin
        definition = ModelDefinition(;
            components=(
                P=Plankton(;
                    states=(:carbon, :chlorophyll), reference_state=:carbon
                ),
                DIC=Pool(:carbon),
            ),
            processes=(
                growth=Growth(;
                    plankton=:P,
                    bindings=(maximum_rate=:maximum_growth_rate,),
                    factors=(
                        nutrients=NutrientResponse(
                            Monod(); resource=:DIC,
                            bindings=(half_saturation=:half_saturation,),
                        ),
                    ),
                ),
            ),
            parameters=(
                maximum_growth_rate=Parameter(0.5),
                half_saturation=Parameter(1.0),
            ),
        )
        model = Agate.Construction.construct(definition)
        _test_tendencies(
            model,
            (P_carbon=2.0, P_chlorophyll=1.0, DIC=1.0),
            (P_carbon=0.5, P_chlorophyll=0.0, DIC=-0.5),
        )
    end

    @testset "non-quota growth rejects additional elemental states" begin
        definition = ModelDefinition(;
            components=(
                P=Plankton(; states=(:carbon, :nitrogen), reference_state=:carbon),
                DIC=Pool(:carbon),
            ),
            processes=(
                growth=Growth(;
                    plankton=:P,
                    factors=(nutrients=NutrientResponse(Monod(); resource=:DIC),),
                ),
            ),
        )
        message = canonicalization_error_message(definition)
        @test occursin("additional elemental states (:nitrogen,)", message)
        @test occursin("NutrientUptake", message)
    end

    @testset "mortality transfers every prognostic state" begin
        definition = ModelDefinition(;
            components=(
                P=Plankton(;
                    states=(:carbon, :nitrogen, :chlorophyll), reference_state=:carbon
                ),
                DOC=Pool(:carbon),
                DON=Pool(:nitrogen),
            ),
            processes=(
                mortality=Mortality(
                    LinearMortality();
                    plankton=:P,
                    bindings=(rate=:mortality_rate,),
                    products=Products((organic=(carbon=:DOC, nitrogen=:DON),)),
                ),
            ),
            parameters=(mortality_rate=Parameter(0.1),),
        )
        model = Agate.Construction.construct(definition)
        state = (
            P_carbon=10.0,
            P_nitrogen=2.0,
            P_chlorophyll=1.0,
            DOC=0.0,
            DON=0.0,
        )
        expected = (
            P_carbon=-1.0,
            P_nitrogen=-0.2,
            P_chlorophyll=-0.1,
            DOC=1.0,
            DON=0.2,
        )
        _test_tendencies(model, state, expected)
    end

    @testset "grazing transfers elemental states and removes physiological states" begin
        definition = ModelDefinition(;
            components=(
                P=Plankton(;
                    states=(:carbon, :nitrogen, :chlorophyll), reference_state=:carbon
                ),
                Z=Plankton(; states=(:carbon, :nitrogen), reference_state=:carbon),
                DOC=Pool(:carbon),
                DON=Pool(:nitrogen),
            ),
            processes=(
                grazing=Consumption(
                    PreferentialGrazing();
                    consumers=:Z,
                    resources=:P,
                    bindings=(
                        maximum_rate=:maximum_rate,
                        half_saturation=:half_saturation,
                        palatability=:palatability,
                        assimilation=:assimilation,
                    ),
                    unassimilated_products=Products((
                        organic=(carbon=:DOC, nitrogen=:DON),
                    )),
                ),
            ),
            parameters=(
                maximum_rate=Parameter(1.0),
                half_saturation=Parameter(10.0),
                palatability=Parameter(1.0),
                assimilation=Parameter(0.5),
            ),
        )
        model = Agate.Construction.construct(definition)
        state = (
            P_carbon=10.0,
            P_nitrogen=2.0,
            P_chlorophyll=1.0,
            Z_carbon=2.0,
            Z_nitrogen=0.4,
            DOC=0.0,
            DON=0.0,
        )
        expected = (
            P_carbon=-1.0,
            P_nitrogen=-0.2,
            P_chlorophyll=-0.1,
            Z_carbon=0.5,
            Z_nitrogen=0.1,
            DOC=0.5,
            DON=0.1,
        )
        _test_tendencies(model, state, expected)
    end
end
