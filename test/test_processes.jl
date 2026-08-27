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
    parameter_bindings,
    driver_identities,
    formulation,
    normalize_model,
    participants,
    PreferentialGrazing

import Agate.Processes: factor_inputs

struct ExternalTestFormulation <: AbstractFormulation end

struct BindingDependencyDefault end
struct ExternalTestFactor{F<:AbstractFormulation} <: AbstractFactor
    formulation::F
end

struct MultiDriverTestFactor{F<:AbstractFormulation} <: AbstractFactor
    formulation::F
end
factor_inputs(::MultiDriverTestFactor) = (
    Agate.Processes.FactorDriver(:wind),
    Agate.Processes.FactorDriver(:temperature),
)

function normalization_error_message(definition)
    err = try
        normalize_model(definition)
        nothing
    catch caught
        caught
    end
    @test err isa ArgumentError
    return err isa Exception ? sprint(showerror, err) : ""
end

@testset "Process authoring and normalization" begin

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

    external_factor = ExternalTestFactor(ExternalTestFormulation())
    external_model = normalize_model(ModelDefinition(;
        components=(P=Population(:nitrogen), N=Pool(:nitrogen)),
        processes=(growth=Growth(;
            populations=:P,
            factors=(light=light, nutrients=response, external=external_factor),
        ),),
    ))
    @test all(binding -> binding.path != (:factors, :external), parameter_bindings(external_model))

    grazing = Consumption(
        PreferentialGrazing();
        consumers=:Z,
        resources=(:P, :B),
        unassimilated_products=:D,
    )
    @test participants(grazing) == (consumer=(:Z,), resource=(:P, :B))

    shared_driver_model = normalize_model(
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

    multi_driver_model = normalize_model(ModelDefinition(;
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
    @test_throws ArgumentError normalize_model(invalid_growth)

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
    @test_throws ArgumentError normalize_model(redundant_growth_source)

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
    @test_throws ArgumentError normalize_model(wrong_currency)

    # Invalid built-in formulation combinations are rejected by their concrete objects.
    @test_throws MethodError Light(Monod(), :PAR, NamedTuple())
    @test_throws MethodError Mortality(Monod(), (:P,), nothing, NamedTuple())
end

@testset "Normalization owns authored structure" begin
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
            "Growth routing",
            one_process(:growth, Growth(; populations=:P, factors=(light=light,)),
                        (P=single, N=Pool(:nitrogen))),
            ("process :growth", "exactly one NutrientResponse or Nutrients"),
        ),
        (
            "Growth maximum-rate owner",
            one_process(:growth, Growth(; populations=:P, factors=(nutrients=nutrient,)),
                        (P=single, N=Pool(:nitrogen))),
            ("process :growth", "exactly one factor that owns the maximum_rate slot"),
        ),
        (
            "Growth ambiguous maximum-rate owner",
            one_process(
                :growth,
                Growth(; populations=:P, factors=(
                    first_light=light,
                    nutrients=nutrient,
                    second_light=Light(Smith(); driver=:PAR),
                )),
                (P=single, N=Pool(:nitrogen)),
            ),
            ("process :growth", "exactly one factor that owns the maximum_rate slot"),
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
            message = normalization_error_message(definition)
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
            message = normalization_error_message(ModelDefinition(;
                components=NamedTuple(), processes=NamedTuple(), parameters
            ))
            @test occursin(fragment, message)
        end
    end
end

@testset "NiPiZD parameter definition validation" begin
    family = Agate.Models.NiPiZD.NiPiZDFamily()
    for parameters in ((), (invalid=1,))
        @test_throws ArgumentError normalize_model(ModelDefinition(;
            components=default_components(family),
            processes=default_processes(family),
            parameters,
        ))
    end
end

@testset "Inline parameter binding resolution" begin
    components = (
        P=Population(:nitrogen; size_structure=[1.0]),
        Z=Population(:nitrogen; size_structure=[10.0]),
        D=Pool(:nitrogen),
        E=Pool(:nitrogen),
    )

    shared = Mortality(
        Agate.Processes.LinearMortality();
        populations=(:P, :Z),
        bindings=(rate=:shared_mortality,),
    )
    remineralization = Agate.Processes.Remineralization(
        Agate.Processes.LinearRemineralization();
        sources=(:D,),
        destination=:D,
        bindings=(rate=(D=:remineralization_rate,),),
    )
    normalized = normalize_model(ModelDefinition(;
        components,
        processes=(mortality=shared, remineralization=remineralization),
        parameters=(
            shared_mortality=Parameter(ConstantDefault(0.1)),
            remineralization_rate=Parameter(ConstantDefault(0.2)),
        ),
    ))

    mortality_bindings = filter(
        binding -> binding.process === :mortality,
        normalized.parameter_bindings,
    )
    @test length(mortality_bindings) == 2
    @test all(binding -> binding.parameter === :shared_mortality, mortality_bindings)
    @test Set(
        (binding.qualifier.axis, binding.qualifier.value) for binding in mortality_bindings
    ) == Set(((:population, :P), (:population, :Z)))
    remineralization_binding = only(filter(
        binding -> binding.process === :remineralization,
        normalized.parameter_bindings,
    ))
    @test (
        remineralization_binding.parameter,
        remineralization_binding.axes,
        (remineralization_binding.qualifier.axis, remineralization_binding.qualifier.value),
    ) == (:remineralization_rate, (), (:source, :D))

    @test_throws MethodError Parameter(ConstantDefault(0.1); axes=:plankton)

    meta_bound_error = normalization_error_message(ModelDefinition(;
        components=(P=components.P, Z=components.Z),
        processes=(mortality=shared,),
        parameters=(shared_mortality=MetaParameter(ConstantDefault(0.1); axes=:plankton),),
    ))
    @test occursin("construction-only and cannot bind to a process slot", meta_bound_error)

    incompatible_axes_error = normalization_error_message(ModelDefinition(;
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
    @test_throws ArgumentError normalize_model(accidental)

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
    @test_throws ArgumentError normalize_model(missing_qualifier)

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
    @test_throws ArgumentError normalize_model(unknown)

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
    @test_throws ArgumentError normalize_model(unknown_zero_slot)

    single_mortality = Mortality(
        Agate.Processes.LinearMortality();
        populations=:P,
        bindings=(rate=:mortality_rate,),
    )
    dependency_only = normalize_model(ModelDefinition(;
        components=(P=components.P,),
        processes=(mortality=single_mortality,),
        parameters=(
            mortality_rate=Parameter(
                DerivedDefault(BindingDependencyDefault(); deps=(:helper,))
            ),
            helper=MetaParameter(ConstantDefault(1.0); axes=:plankton),
        ),
    ))
    @test dependency_only.parameters.helper.axes === :plankton
    @test dependency_only.parameters.mortality_rate.default.deps == (:helper,)

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
        message = normalization_error_message(ModelDefinition(;
            components=(P=components.P,),
            processes=(mortality=single_mortality,),
            parameters,
        ))
        @test occursin(fragment, message)
    end
end
