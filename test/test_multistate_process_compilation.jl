using Test

using Agate.Components: Plankton, Pool
using Agate.Parameters: Parameter
using Agate.Processes:
    Consumption, Growth, LinearMortality, Monod, Mortality, ModelDefinition,
    NutrientResponse, PreferentialGrazing, Products, QuadraticMortality

@testset "Multi-state process semantics" begin
    @testset "growth factors are optional" begin
        definition = ModelDefinition(;
            components=(
                P=Plankton(; states=(carbon=:carbon,), reference_state=:carbon),
                DIC=Pool(:carbon),
            ),
            processes=(
                growth=Growth(;
                    plankton=:P,
                    reference_resource=:DIC,
                    bindings=(maximum_rate=:maximum_growth_rate,),
                ),
            ),
            parameters=(maximum_growth_rate=Parameter(0.5),),
        )
        model = Agate.Construction.construct(definition)
        test_tendencies(model, (P=2.0, DIC=10.0), (P=1.0, DIC=-1.0))
    end

    @testset "growth permits non-elemental physiological states" begin
        definition = ModelDefinition(;
            components=(
                P=Plankton(;
                    states=(carbon=:carbon, chlorophyll=nothing), reference_state=:carbon
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

    @testset "component composition drives elemental growth" begin
        definition = ModelDefinition(;
            components=(
                P=Plankton(; states=(carbon=:carbon,), reference_state=:carbon,
                    size_structure=[1.0, 2.0], composition=(nitrogen=:nitrogen_to_carbon,)),
                DIC=Pool(:carbon),
                DIN=Pool(:nitrogen),
            ),
            processes=(
                growth=Growth(;
                    plankton=:P,
                    reference_resource=:DIC,
                    additional_resources=(nitrogen=:DIN,),
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
                nitrogen_to_carbon=Parameter([0.2, 0.4]),
            ),
        )
        model = Agate.Construction.construct(definition)
        test_tendencies(
            model,
            (P_1=2.0, P_2=3.0, DIC=10.0, DIN=10.0),
            (P_1=1.0, P_2=1.5, DIC=-2.5, DIN=-0.8),
        )
        @test_throws ArgumentError Agate.Construction.construct(
            definition; parameter_overrides=(nitrogen_to_carbon=[0.2, -0.1],)
        )
    end

    @testset "growth products close every implicit element" begin
        components = (
            P=Plankton(; states=(carbon=:carbon,), reference_state=:carbon,
                composition=(nitrogen=:nitrogen_to_carbon,)),
            DIC=Pool(:carbon), DIN=Pool(:nitrogen), DOC=Pool(:carbon), DON=Pool(:nitrogen),
        )
        growth(products) = Growth(;
            plankton=:P, reference_resource=:DIC, additional_resources=(nitrogen=:DIN,),
            products,
            bindings=(maximum_rate=:maximum_growth_rate, product_fraction=:exudation_fraction),
        )
        parameters = (
            maximum_growth_rate=Parameter(0.5), nitrogen_to_carbon=Parameter(0.2),
            exudation_fraction=Parameter(0.2),
        )

        products = Products((exudate=(carbon=:DOC, nitrogen=:DON),))
        model = Agate.Construction.construct(ModelDefinition(;
            components, processes=(growth=growth(products),), parameters,
        ))
        test_tendencies(
            model, (P=2.0, DIC=10.0, DIN=10.0, DOC=0.0, DON=0.0),
            (P=0.8, DIC=-1.0, DIN=-0.2, DOC=0.2, DON=0.04),
        )

        for invalid in (
            Products(:DOC),
            Products((exudate=(carbon=:DOC,),)),
        )
            @test_throws ArgumentError Agate.Construction.construct(ModelDefinition(;
                components, processes=(growth=growth(invalid),), parameters,
            ))
        end
    end

    @testset "growth products require sourced elements" begin
        products = Products((exudate=(nitrogen=:DON, carbon=:DOC),))
        parameters = (maximum_growth_rate=Parameter(0.5), product_fraction=Parameter(0.25))
        growth = Growth(; plankton=:P, reference_resource=:DIN, products,
            bindings=(maximum_rate=:maximum_growth_rate, product_fraction=:product_fraction))
        for states in ((nitrogen=:nitrogen,), (nitrogen=:nitrogen, carbon=:carbon))
            components = (
                P=Plankton(; states, reference_state=:nitrogen),
                DIN=Pool(:nitrogen), DON=Pool(:nitrogen), DOC=Pool(:carbon),
            )
            message = canonicalization_error_message(ModelDefinition(;
                components, processes=(; growth), parameters,
            ))
            @test occursin("growth products product elements", message)
        end

        components = (
            P=Plankton(; states=(nitrogen=:nitrogen,), reference_state=:nitrogen),
            DIN=Pool(:nitrogen), DON=Pool(:nitrogen),
        )
        definition = ModelDefinition(;
            components,
            processes=(growth=Growth(; plankton=:P, reference_resource=:DIN,
                products=Products((exudate=:DON,)),
                bindings=(maximum_rate=:maximum_growth_rate, product_fraction=:product_fraction)),),
            parameters,
        )
        model = Agate.Construction.construct(definition)
        test_tendencies(model, (P=2.0, DIN=10.0, DON=0.0),
            (P=0.75, DIN=-1.0, DON=0.25))
    end

    @testset "heterogeneous prey preserve shared grazing and elemental budgets" begin
        components = (
            P1=Plankton(; states=(carbon=:carbon,), reference_state=:carbon,
                composition=(nitrogen=:p1_n2c,)),
            P2=Plankton(; states=(carbon=:carbon,), reference_state=:carbon,
                composition=(nitrogen=:p2_n2c,)),
            Z=Plankton(; states=(carbon=:carbon, nitrogen=:nitrogen), reference_state=:carbon),
            DOC=Pool(:carbon), DON=Pool(:nitrogen),
        )
        processes = (grazing=Consumption(PreferentialGrazing();
            consumers=:Z, resources=(:P1, :P2),
            unassimilated_products=Products((organic=(carbon=:DOC, nitrogen=:DON),)),
        ),)
        parameters = (
            maximum_rate=Parameter(1.0), half_saturation=Parameter(1.0),
            palatability=Parameter(1.0), assimilation=Parameter(0.5),
            p1_n2c=Parameter(0.2), p2_n2c=Parameter(0.4),
        )
        model = Agate.Construction.construct(ModelDefinition(; components, processes, parameters))
        state = (P1=1.0, P2=1.0, Z_carbon=1.0, Z_nitrogen=0.1, DOC=0.0, DON=0.0)
        names = Agate.Introspection.tracer_names(model)
        args = (0, 0, 0, 0, (getproperty(state, name) for name in names)...)
        t = model_tendencies(model, args)
        @test t.P1 < 0 && t.P2 < 0
        @test isapprox(t.P1 + t.P2 + t.Z_carbon + t.DOC, 0; atol=1e-12)
        @test isapprox(0.2 * t.P1 + 0.4 * t.P2 + t.Z_nitrogen + t.DON, 0; atol=1e-12)
    end

    @testset "mortality routes mixed explicit and implicit elemental inventories" begin
        components = (
            P=Plankton(; states=(carbon=:carbon, nitrogen=:nitrogen),
                reference_state=:carbon, composition=(phosphorus=:p2c,)),
            DOC=Pool(:carbon), DON=Pool(:nitrogen), DOP=Pool(:phosphorus),
        )
        processes = (mortality=Mortality(LinearMortality(); plankton=:P,
            bindings=(rate=:mortality_rate,),
            products=Products((organic=(carbon=:DOC, nitrogen=:DON, phosphorus=:DOP),))),)
        parameters = (mortality_rate=Parameter(0.1), p2c=Parameter(0.1))
        model = Agate.Construction.construct(ModelDefinition(; components, processes, parameters))
        test_tendencies(model,
            (P_carbon=10.0, P_nitrogen=2.0, DOC=0.0, DON=0.0, DOP=0.0),
            (P_carbon=-1.0, P_nitrogen=-0.2, DOC=1.0, DON=0.2, DOP=0.1))
    end

    @testset "growth leaves independently prognostic elemental states unchanged" begin
        definition = ModelDefinition(;
            components=(
                P=Plankton(; states=(biomass=:carbon, reserve_n=:nitrogen), reference_state=:biomass),
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
            (DIC=10.0, P_biomass=2.0, P_reserve_n=0.3),
            (DIC=-1.0, P_biomass=1.0, P_reserve_n=0.0),
        )
    end

    @testset "mortality transfers every prognostic state" begin
        components = (
            P=Plankton(;
                states=(carbon=:carbon, nitrogen=:nitrogen, chlorophyll=nothing), reference_state=:carbon
            ),
            DOC=Pool(:carbon),
            DON=Pool(:nitrogen),
        )
        state = (
            P_carbon=10.0, P_nitrogen=2.0, P_chlorophyll=1.0, DOC=0.0, DON=0.0,
        )
        expected = (
            P_carbon=-1.0, P_nitrogen=-0.2, P_chlorophyll=-0.1, DOC=1.0, DON=0.2,
        )

        # Match reference-state loss intensity (0.1). The quadratic case then
        # verifies that every state uses carbon intensity, not its own inventory squared.
        for (formulation, mortality_rate) in (
            (LinearMortality(), 0.1),
            (QuadraticMortality(), 0.01),
        )
            mortality = Mortality(
                formulation;
                plankton=:P,
                bindings=(rate=:mortality_rate,),
                products=Products((organic=(carbon=:DOC, nitrogen=:DON),)),
            )
            definition = ModelDefinition(;
                components,
                processes=(; mortality),
                parameters=(mortality_rate=Parameter(mortality_rate),),
            )
            model = Agate.Construction.construct(definition)
            test_tendencies(model, state, expected)
        end
    end

    @testset "grazing transfers elemental states and removes physiological states" begin
        definition = ModelDefinition(;
            components=(
                P=Plankton(;
                    states=(carbon=:carbon, nitrogen=:nitrogen, chlorophyll=nothing), reference_state=:carbon
                ),
                Z=Plankton(; states=(carbon=:carbon, nitrogen=:nitrogen), reference_state=:carbon),
                DOC=Pool(:carbon),
                DON=Pool(:nitrogen),
            ),
            processes=(
                grazing=Consumption(
                    PreferentialGrazing(; switching_exponent=2);
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

@testset "Grouped processes route heterogeneous elemental inventories" begin
    components = (
        P1=Plankton(; states=(carbon=:carbon,), reference_state=:carbon),
        P2=Plankton(; states=(carbon=:carbon,), reference_state=:carbon,
            composition=(nitrogen=:p2_n2c,)),
        Z=Plankton(; states=(carbon=:carbon, nitrogen=:nitrogen), reference_state=:carbon),
        DIC=Pool(:carbon), DIN=Pool(:nitrogen), DOC=Pool(:carbon), DON=Pool(:nitrogen),
    )
    products = Products((organic=(carbon=:DOC, nitrogen=:DON),))
    build(process, parameters) = Agate.Construction.construct(ModelDefinition(;
        components, processes=(transfer=process,),
        parameters=merge((p2_n2c=Parameter(0.2),), parameters),
    ))
    state = (P1=2.0, P2=3.0, Z_carbon=1.0, Z_nitrogen=0.1,
        DIC=10.0, DIN=10.0, DOC=0.0, DON=0.0)

    growth = Growth(; plankton=(:P1, :P2), reference_resource=:DIC,
        additional_resources=(nitrogen=:DIN,), products,
        bindings=(maximum_rate=:growth_rate, product_fraction=:product_fraction))
    test_tendencies(build(growth, (growth_rate=Parameter(0.5), product_fraction=Parameter(0.2))),
        state, (P1=0.8, P2=1.2, DIC=-2.5, DIN=-0.3, DOC=0.5, DON=0.06))

    mortality = Mortality(LinearMortality(); plankton=(:P1, :P2), products,
        bindings=(rate=:mortality_rate,))
    test_tendencies(build(mortality, (mortality_rate=Parameter(0.1),)), state,
        (P1=-0.2, P2=-0.3, DOC=0.5, DON=0.06))

    grazing = Consumption(PreferentialGrazing(); consumers=:Z, resources=(:P1, :P2),
        unassimilated_products=products)
    test_tendencies(build(grazing, (
        maximum_rate=Parameter(1.0), half_saturation=Parameter(1.0),
        palatability=Parameter(1.0), assimilation=Parameter(0.5),
    )), state, (P1=-1 / 3, P2=-1 / 2, Z_carbon=5 / 12, Z_nitrogen=0.05,
        DOC=5 / 12, DON=0.05))

    invalid = ModelDefinition(; components,
        processes=(transfer=Mortality(LinearMortality(); plankton=(:P1, :P2),
            products=Products((organic=(carbon=:DOC,),))),),
        parameters=(p2_n2c=Parameter(0.2), mortality_rate=Parameter(0.1)),
    )
    @test occursin("product elements", canonicalization_error_message(invalid))
end
