using Test

using Agate.Configuration: Plankton, Pool
using Agate.Parameters: Parameter
using Agate.Processes:
    Consumption, FixedStoichiometry, Growth, LinearMortality, Monod, Mortality, ModelDefinition,
    NutrientResponse, PreferentialGrazing, Products

@testset "Multi-state process semantics" begin
    @testset "growth permits non-elemental physiological states" begin
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
                    reference_resource=:DIC,
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
        test_tendencies(
            model,
            (P_carbon=2.0, P_chlorophyll=1.0, DIC=1.0),
            (P_carbon=0.5, P_chlorophyll=0.0, DIC=-0.5),
        )
    end

    @testset "fixed-stoichiometry transfer is process-owned" begin
        definition = ModelDefinition(;
            components=(
                P=Plankton(; states=:carbon, reference_state=:carbon),
                DIC=Pool(:carbon),
                DIN=Pool(:nitrogen),
            ),
            processes=(
                growth=Growth(;
                    plankton=:P,
                    reference_resource=:DIC,
                    additional_resources=(nitrogen=:DIN,),
                    stoichiometry=FixedStoichiometry(;
                        reference_element=:carbon,
                        bindings=(ratio=(nitrogen=:nitrogen_to_carbon,),),
                    ),
                    bindings=(maximum_rate=:maximum_growth_rate,),
                    factors=(
                        carbon=NutrientResponse(
                            Monod(); resource=:DIC, bindings=(half_saturation=:half_saturation,)
                        ),
                    ),
                ),
            ),
            parameters=(
                maximum_growth_rate=Parameter(0.5),
                half_saturation=Parameter(0.0),
                nitrogen_to_carbon=Parameter(0.2),
            ),
        )
        model = Agate.Construction.construct(definition)
        test_tendencies(
            model,
            (P=2.0, DIC=10.0, DIN=10.0),
            (P=1.0, DIC=-1.0, DIN=-0.2),
        )
    end

    @testset "growth leaves independently prognostic elemental states unchanged" begin
        definition = ModelDefinition(;
            components=(
                P=Plankton(; states=(:carbon, :nitrogen), reference_state=:carbon),
                DIC=Pool(:carbon),
            ),
            processes=(
                growth=Growth(;
                    plankton=:P,
                    reference_resource=:DIC,
                    bindings=(maximum_rate=:maximum_growth_rate,),
                    factors=(
                        carbon=NutrientResponse(
                            Monod(); resource=:DIC,
                            bindings=(half_saturation=:half_saturation,),
                        ),
                    ),
                ),
            ),
            parameters=(
                maximum_growth_rate=Parameter(0.5),
                half_saturation=Parameter(0.0),
            ),
        )
        model = Agate.Construction.construct(definition)
        test_tendencies(
            model,
            (DIC=10.0, P_carbon=2.0, P_nitrogen=0.3),
            (DIC=-1.0, P_carbon=1.0, P_nitrogen=0.0),
        )
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
        test_tendencies(model, state, expected)
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
        test_tendencies(model, state, expected)
    end
end
