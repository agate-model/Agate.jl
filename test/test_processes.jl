using Test
using Agate.Configuration: Population, Pool
using Agate.ModelFamilies: default_components, default_processes
using Agate.Parameters: ConstantDefault, DerivedDefault, MetaParameter, Parameter
using Agate.Processes:
    AbstractFormulation,
    AbstractFactor,
    Smith,
    Geider,
    Monod,
    Liebig,
    MultiplicativeFactors,
    Growth,
    Light,
    Nutrients,
    NutrientResponse,
    FixedStoichiometry,
    Consumption,
    Mortality,
    ModelDefinition,
    driver_identities,
    formulation,
    HeterotrophicConsumption,
    canonicalize_model,
    participants,
    PreferentialGrazing

import Agate.Processes: factor_inputs

struct ExternalTestFormulation <: AbstractFormulation end

struct BindingDependencyDefault end
struct MultiDriverTestFactor{F<:AbstractFormulation} <: AbstractFactor
    formulation::F
end
factor_inputs(::MultiDriverTestFactor) = (
    Agate.Processes.FactorDriver(:wind),
    Agate.Processes.FactorDriver(:temperature),
)

@testset "Process authoring and canonicalization" begin

    light = Light(Smith(); driver=:PAR)
    response = NutrientResponse(Monod(); resource=:N)
    growth = Growth(;
        populations=:P,
        factors=(light=light, nutrients=response),
    )
    @test formulation(light) isa Smith
    @test formulation(response) isa Monod
    @test formulation(growth) isa MultiplicativeFactors

    @test participants(growth) == (population=(:P,),)

    grazing = Consumption(
        PreferentialGrazing();
        consumers=:Z,
        resources=(:P, :B),
        unassimilated_products=:D,
    )
    @test participants(grazing) == (consumer=(:Z,), resource=(:P, :B))

    shared_driver_model = canonicalize_model(
        ModelDefinition(;
            components=(
                P=Population(:nitrogen), Z=Population(:nitrogen), N=Pool(:nitrogen)
            ),
            processes=(
                growth_Z=Growth(;
                    populations=:Z,
                    factors=(
                        light=Light(Smith(); driver=:PAR),
                        nutrients=NutrientResponse(Monod(); resource=:N),
                    ),
                ),
                growth_P=Growth(;
                    populations=:P,
                    factors=(
                        light=Light(Smith(); driver=:PAR),
                        nutrients=NutrientResponse(Monod(); resource=:N),
                    ),
                ),
            ),
        ),
    )
    @test driver_identities(shared_driver_model) == (:PAR,)

    multi_driver_model = canonicalize_model(ModelDefinition(;
        components=(P=Population(:nitrogen), N=Pool(:nitrogen)),
        processes=(growth=Growth(;
            populations=:P,
            factors=(
                light=Light(Smith(); driver=:PAR),
                nutrients=NutrientResponse(Monod(); resource=:N),
                environment=MultiDriverTestFactor(ExternalTestFormulation()),
            ),
        ),),
    ))
    @test driver_identities(multi_driver_model) == (:PAR, :temperature, :wind)

    invalid_growth = ModelDefinition(;
        components=(P=Population(:nitrogen),),
        processes=(growth=Growth(;
            populations=:P,
            factors=(nutrients=NutrientResponse(Monod(); resource=:missing),),
        ),),
    )
    @test_throws ArgumentError canonicalize_model(invalid_growth)

    @test_throws MethodError Light(:smith; driver=:PAR)
    @test_throws ArgumentError Growth(; populations=:P, factors=NamedTuple())
    @test_throws ArgumentError Growth(;
        populations=(), factors=(light=Light(Smith(); driver=:PAR),)
    )
    @test_throws ArgumentError Consumption(PreferentialGrazing(); consumers=(), resources=:P)
    @test_throws ArgumentError Consumption(PreferentialGrazing(); consumers=:Z, resources=(:P, 1))

    redundant_growth_source = ModelDefinition(;
        components=(
            P=Population(:nitrogen),
            N=Pool(:nitrogen),
            D=Pool(:nitrogen),
        ),
        processes=(growth=Growth(;
            populations=:P,
            source=:D,
            factors=(
                light=Light(Smith(); driver=:PAR),
                nutrients=NutrientResponse(Monod(); resource=:N),
            ),
        ),),
    )
    @test_throws ArgumentError canonicalize_model(redundant_growth_source)

    wrong_currency = ModelDefinition(;
        components=(
            P=Population(:carbon),
            DIC=Pool(:carbon),
            DIN=Pool(:phosphorus),
        ),
        processes=(growth=Growth(;
            populations=:P,
            source=:DIC,
            factors=(
                light=Light(Geider(); driver=:PAR),
                nutrients=Nutrients(
                    Liebig();
                    responses=(nitrogen=NutrientResponse(Monod(); resource=:DIN),),
                ),
            ),
            stoichiometry=FixedStoichiometry(; reference=:carbon),
        ),),
    )
    @test_throws ArgumentError canonicalize_model(wrong_currency)

    # Invalid built-in formulation combinations are rejected by their concrete objects.
    @test_throws MethodError Light(Monod(), :PAR, NamedTuple())
    @test_throws MethodError Mortality(Monod(), (:P,), nothing, NamedTuple())
end

@testset "Canonicalization owns authored structure" begin
    single = Population(:nitrogen)
    multi = Population(; states=(carbon=:carbon, nitrogen=:nitrogen))
    light = Light(Smith(); driver=:PAR)
    nutrient = NutrientResponse(Monod(); resource=:N)
    multi_nutrient = Nutrients(
        Liebig(); responses=(nitrogen=NutrientResponse(Monod(); resource=:N),)
    )
    one_process(id, process, components) = ModelDefinition(;
        components, processes=NamedTuple{(id,)}((process,))
    )

    cases = (
        (
            "Growth nutrient factor",
            one_process(:growth, Growth(; populations=:P, factors=(light=light,)),
                        (P=single, N=Pool(:nitrogen))),
            ("process :growth", "exactly one NutrientResponse or Nutrients"),
        ),
        (
            "multi-resource Growth source",
            one_process(:growth, Growth(;
                populations=:P, factors=(light=light, nutrients=multi_nutrient)
            ), (P=Population(:carbon), N=Pool(:nitrogen))),
            ("process :growth", "requires a source component"),
        ),
        (
            "multi-resource Growth stoichiometry",
            one_process(:growth, Growth(;
                populations=:P, source=:DIC, factors=(light=light, nutrients=multi_nutrient)
            ), (P=Population(:carbon), DIC=Pool(:carbon), N=Pool(:nitrogen))),
            ("process :growth", "requires FixedStoichiometry"),
        ),
        (
            "structured nutrient resource",
            one_process(:growth, Growth(;
                populations=:P, factors=(light=light, nutrients=nutrient)
            ), (P=single, N=Pool(:nitrogen; size_structure=[1.0]))),
            ("process :growth", "nutrient factor resource", "scalar Pool"),
        ),
        (
            "structured factor resource",
            one_process(
                :consume,
                Consumption(
                    HeterotrophicConsumption();
                    consumers=:B,
                    resources=:POM,
                    factors=(limit=NutrientResponse(Monod(); resource=:POM),),
                ),
                (B=Population(:nitrogen), POM=Pool(:nitrogen; size_structure=[1.0])),
            ),
            ("process :consume", "component :POM", "scalar component"),
        ),
        (
            "multi-state Growth",
            one_process(:growth, Growth(;
                populations=:P, factors=(light=light, nutrients=nutrient)
            ), (P=multi, N=Pool(:carbon))),
            ("process :growth", "explicit state selection"),
        ),
        (
            "multi-state Consumption",
            one_process(:consume, Consumption(
                PreferentialGrazing(); consumers=:Z, resources=:P
            ), (Z=multi, P=Population(:carbon))),
            ("process :consume", "explicit state selection"),
        ),
        (
            "multi-state Mortality",
            one_process(:mortality, Mortality(
                Agate.Processes.LinearMortality(); populations=:P
            ), (P=multi,)),
            ("process :mortality", "explicit state selection"),
        ),
        (
            "Mortality participant type",
            one_process(:mortality, Mortality(
                Agate.Processes.LinearMortality(); populations=:D
            ), (D=Pool(:nitrogen),)),
            ("process :mortality", "must be a Population"),
        ),
        (
            "Remineralization currency",
            one_process(:remineralization, Agate.Processes.Remineralization(
                Agate.Processes.LinearRemineralization(); sources=:D, destination=:N
            ), (D=Pool(:nitrogen), N=Pool(:phosphorus))),
            ("process :remineralization", "remineralization source :D", "currency"),
        ),
        (
            "structured Remineralization source",
            one_process(:remineralization, Agate.Processes.Remineralization(
                Agate.Processes.LinearRemineralization(); sources=:D, destination=:N
            ), (D=Pool(:nitrogen; size_structure=[1.0]), N=Pool(:nitrogen))),
            ("process :remineralization", "source :D", "scalar Pool"),
        ),
        (
            "structured Remineralization destination",
            one_process(:remineralization, Agate.Processes.Remineralization(
                Agate.Processes.LinearRemineralization(); sources=:D, destination=:N
            ), (D=Pool(:nitrogen), N=Pool(:nitrogen; size_structure=[1.0]))),
            ("process :remineralization", "destination :N", "scalar Pool"),
        ),
        (
            "structured product target",
            one_process(:mortality, Mortality(
                Agate.Processes.LinearMortality(); populations=:P, products=:D
            ), (P=single, D=Pool(:nitrogen; size_structure=[1.0]))),
            ("process :mortality", "product", "scalar Pool"),
        ),
    )

    for (label, definition, fragments) in cases
        @testset "$label" begin
            message = canonicalization_error_message(definition)
            @test all(fragment -> occursin(fragment, message), fragments)
        end
    end

    parameter_cases = (
        ("reserved parameter", (x=Parameter(ConstantDefault(1.0)),), "reserved parameter name :x"),
        (
            "unknown dependency",
            (a=Parameter(DerivedDefault(BindingDependencyDefault(); deps=(:missing,))),),
            "depends on undeclared parameter :missing",
        ),
        (
            "derived dependency",
            (
                a=Parameter(DerivedDefault(BindingDependencyDefault(); deps=(:b,))),
                b=Parameter(DerivedDefault(BindingDependencyDefault(); deps=(:base,))),
                base=Parameter(ConstantDefault(1.0)),
            ),
            "depends on derived parameter :b",
        ),
    )
    for (label, parameters, fragment) in parameter_cases
        @testset "$label" begin
            message = canonicalization_error_message(ModelDefinition(;
                components=NamedTuple(), processes=NamedTuple(), parameters
            ))
            @test occursin(fragment, message)
        end
    end
end

@testset "NiPiZD parameter definition validation" begin
    family = Agate.Models.NiPiZD.NiPiZDFamily()
    for parameters in ((), (invalid=1,))
        @test_throws ArgumentError canonicalize_model(ModelDefinition(;
            components=default_components(family),
            processes=default_processes(family),
            parameters,
        ))
    end
end

@testset "Parameter binding behavior and validation" begin
    components = (
        P=Population(:nitrogen; size_structure=[1.0]),
        Z=Population(:nitrogen; size_structure=[10.0]),
        D=Pool(:nitrogen),
        E=Pool(:nitrogen),
        R=Pool(:nitrogen),
    )

    qualified = Agate.Processes.Remineralization(
        Agate.Processes.LinearRemineralization();
        sources=(:D, :E),
        destination=:R,
        bindings=(rate=(D=:D_rate, E=:E_rate),),
    )
    qualified_model = Agate.Construction.construct(ModelDefinition(;
        components=(D=components.D, E=components.E, R=components.R),
        processes=(remineralization=qualified,),
        parameters=(D_rate=Parameter(0.1), E_rate=Parameter(0.2)),
    ))
    tracer_order = Agate.Introspection.tracer_names(qualified_model)
    state = (D=2.0, E=3.0, R=0.0)
    args = (0.0, 0.0, 0.0, 0.0, Tuple(getproperty(state, tracer) for tracer in tracer_order)...)
    @test qualified_model(Val(:D), args...) ≈ -0.2
    @test qualified_model(Val(:E), args...) ≈ -0.6
    @test qualified_model(Val(:R), args...) ≈ 0.8

    shared = Mortality(
        Agate.Processes.LinearMortality();
        populations=(:P, :Z),
        bindings=(rate=:shared_mortality,),
    )
    @test_throws MethodError Parameter(ConstantDefault(0.1); axes=:plankton)
    @test_throws ArgumentError MetaParameter(ConstantDefault(0.1); axes=:consumer)

    meta_bound_error = canonicalization_error_message(ModelDefinition(;
        components=(P=components.P, Z=components.Z),
        processes=(mortality=shared,),
        parameters=(shared_mortality=MetaParameter(ConstantDefault(0.1); axes=:plankton),),
    ))
    @test occursin("construction-only and cannot bind to a process slot", meta_bound_error)

    incompatible_axes_error = canonicalization_error_message(ModelDefinition(;
        components=(P=components.P, Z=components.Z, D=components.D),
        processes=(
            mortality=shared,
            grazing=Consumption(
                PreferentialGrazing();
                consumers=:Z,
                resources=:P,
                bindings=(maximum_rate=:shared_mortality,),
                unassimilated_products=:D,
            ),
        ),
        parameters=(
            shared_mortality=Parameter(0.1),
            half_saturation=Parameter(0.1),
            palatability=Parameter(1.0),
            assimilation=Parameter(0.7),
        ),
    ))
    @test occursin("incompatible semantic axes", incompatible_axes_error)

    accidental = ModelDefinition(;
        components=(P=components.P, Z=components.Z),
        processes=(
            mortality=Mortality(
                Agate.Processes.LinearMortality(); populations=(:P, :Z)
            ),
        ),
        parameters=(rate=Parameter(ConstantDefault(0.1)),),
    )
    @test_throws ArgumentError canonicalize_model(accidental)

    missing_qualifier = ModelDefinition(;
        components=(D=components.D, E=components.E),
        processes=(
            remineralization=Agate.Processes.Remineralization(
                Agate.Processes.LinearRemineralization();
                sources=(:D, :E),
                destination=:D,
                bindings=(rate=(D=:remineralization_rate,),),
            ),
        ),
        parameters=(remineralization_rate=Parameter(ConstantDefault(0.2)),),
    )
    @test_throws ArgumentError canonicalize_model(missing_qualifier)

    unknown = ModelDefinition(;
        components=(P=components.P,),
        processes=(
            mortality=Mortality(
                Agate.Processes.LinearMortality();
                populations=:P,
                bindings=(missing=:shared_mortality,),
            ),
        ),
        parameters=(shared_mortality=Parameter(ConstantDefault(0.1)),),
    )
    @test_throws ArgumentError canonicalize_model(unknown)

    unknown_zero_slot = ModelDefinition(;
        components=(P=components.P, D=components.D),
        processes=(
            growth=Growth(;
                populations=:P,
                factors=(
                    light=Light(Smith(); driver=:PAR),
                    nutrients=Nutrients(
                        Liebig();
                        responses=(nitrogen=NutrientResponse(Monod(); resource=:D),),
                        bindings=(missing=:shared_mortality,),
                    ),
                ),
            ),
        ),
        parameters=(K=Parameter(ConstantDefault(0.1)),),
    )
    @test_throws ArgumentError canonicalize_model(unknown_zero_slot)

    single_mortality = Mortality(
        Agate.Processes.LinearMortality();
        populations=:P,
        bindings=(rate=:mortality_rate,),
    )
    lifecycle_cases = (
        (
            (mortality_rate=Parameter(0.1), unused=Parameter(1.0)),
            "use MetaParameter for construction-only",
        ),
        (
            (mortality_rate=Parameter(0.1), helper=MetaParameter(1.0)),
            "must be used by at least one DerivedDefault",
        ),
    )
    for (parameters, fragment) in lifecycle_cases
        message = canonicalization_error_message(ModelDefinition(;
            components=(P=components.P,),
            processes=(mortality=single_mortality,),
            parameters,
        ))
        @test occursin(fragment, message)
    end
end
