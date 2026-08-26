using Test
using ForwardDiff

using Agate.Configuration: Population, Pool, realize_components
using Agate.ModelFamilies: default_components, default_processes
using Agate.Parameters: ConstantDefault, DerivedDefault, Parameter
using Agate.Processes:
    AbstractFormulation,
    AbstractFactor,
    Smith,
    Geider,
    Monod,
    Liebig,
    FrankTNorm,
    Q10,
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
    resolve_parameter_applicability,
    driver_identities,
    formulation,
    factor_value,
    normalize_model,
    participants,
    process_rate,
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
    @test :Light in names(Agate.Processes)
    @test :FixedStoichiometry in names(Agate.Processes)
    @test :Grazing ∉ names(Agate.Processes)
    @test :Smith in names(Agate.Processes)

    light = Light(Smith(); driver=:PAR)
    response = NutrientResponse(Monod(); resource=:N)
    growth = Growth(;
        populations=:P,
        factors=(light=light, nutrients=response),
    )
    @test formulation(light) isa Smith
    @test formulation(response) isa Monod
    @test formulation(growth) isa MultiplicativeFactors

    frank_nutrients = Nutrients(
        FrankTNorm(); responses=(nitrogen=NutrientResponse(Monod(); resource=:N),)
    )
    @test fieldcount(FrankTNorm) == 0
    @test participants(growth) == (population=(:P,),)
    @test :drivers ∉ names(Agate.Processes)

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
            "dependency cycle",
            (
                a=Parameter(DerivedDefault(BindingDependencyDefault(); deps=(:b,))),
                b=Parameter(DerivedDefault(BindingDependencyDefault(); deps=(:a,))),
            ),
            "dependency cycle",
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

@testset "Compositional factor mathematics" begin
    biomass = 0.4
    maximum_rate = 2e-5
    light = 100.0
    alpha = 2e-6

    base = process_rate(MultiplicativeFactors(), biomass, maximum_rate)
    smith = factor_value(Smith(), light, maximum_rate, alpha)
    monod = factor_value(Monod(), 0.5, 0.2)
    @test isapprox(
        base * smith * monod,
        maximum_rate * biomass *
        (alpha * light / sqrt(maximum_rate^2 + (alpha * light)^2)) *
        (0.5 / 0.7);
        rtol=1e-14,
    )

    chlorophyll_to_carbon = 0.02
    geider = factor_value(
        Geider(), light, maximum_rate, alpha, chlorophyll_to_carbon
    )
    @test geider == Agate.Library.Photosynthesis.geider_light_response(
        light, alpha, maximum_rate, chlorophyll_to_carbon
    )
    @test factor_value(Geider(), light, 0.0, alpha, chlorophyll_to_carbon) == 0.0
    limitations = (factor_value(Monod(), 3.0, 0.5), factor_value(Monod(), 0.2, 0.5))
    liebig = factor_value(Liebig(), limitations)
    frank = factor_value(FrankTNorm(), limitations, 25)
    @test frank ≈ Agate.Library.Nutrients.frank_tnorm(limitations; sharpness=25)
    frank_gradient = ForwardDiff.gradient(
        x -> factor_value(FrankTNorm(), (x[1], x[2]), 25), collect(limitations)
    )
    @test all(isfinite, frank_gradient)
    temperature = factor_value(Q10(), 30.0, 2.0, 20.0)
    @test temperature == Agate.Library.Temperature.q10_temperature_factor(30.0, 2.0, 20.0)
    expected_liebig = min(3.0 / 3.5, 0.2 / 0.7)
    expected_geider = 1 - exp(-alpha * chlorophyll_to_carbon * light / maximum_rate)
    @test isapprox(
        base * geider * liebig * temperature,
        maximum_rate * biomass * expected_geider * expected_liebig * 2;
        rtol=1e-14,
    )
end

@testset "NiPiZD normalized process contract" begin
    family = Agate.Models.NiPiZD.NiPiZDFamily()
    definition = ModelDefinition(family)
    normalized = normalize_model(definition)

    @test normalized.components == default_components(family)
    growth_facts = normalized.processes.growth_P.facts
    @test growth_facts.routing.mode === :single_resource
    @test growth_facts.routing.factor === :nutrients
    @test growth_facts.maximum_rate_factor === :light
    @test Tuple((ref.population, ref.state) for ref in growth_facts.population_states) ==
        ((:P, :nitrogen),)

    processes = default_processes(family)
    reversed_names = reverse(Tuple(keys(processes)))
    reversed_processes = NamedTuple{reversed_names}(
        Tuple(getproperty(processes, name) for name in reversed_names)
    )
    reordered = normalize_model(ModelDefinition(;
        components=default_components(family),
        processes=reversed_processes,
        parameters=definition.parameters,
    ))
    binding_key(binding) = (
        binding.process, binding.path, binding.slot,
        isnothing(binding.qualifier) ? nothing : (binding.qualifier.axis, binding.qualifier.value),
        binding.parameter,
    )
    @test map(binding_key, reordered.parameter_bindings) ==
        map(binding_key, normalized.parameter_bindings)

    layout = realize_components(normalized.components)
    applicability = resolve_parameter_applicability(normalized, layout)
    function application(parameter, process, slot; path=())
        return only(filter(applicability) do item
            binding = item.binding
            binding.parameter === parameter && binding.process === process &&
                binding.path == path && binding.slot === slot
        end)
    end
    @test application(
        :maximum_growth_rate, :growth_P, :maximum_rate; path=(:factors, :light)
    ).axis_classes == ((:P_1, :P_2),)
    @test application(:palatability_matrix, :grazing_Z_on_P, :palatability).axis_classes ==
        ((:Z_1, :Z_2), (:P_1, :P_2))

    invalid = ModelDefinition(;
        components=default_components(family),
        processes=(bad=Consumption(PreferentialGrazing(); consumers=:Z, resources=:missing),),
    )
    @test_throws ArgumentError normalize_model(invalid)

    @test_throws ArgumentError normalize_model(
        ModelDefinition(;
            components=default_components(family),
            processes=default_processes(family),
            parameters=(),
        ),
    )
    @test_throws ArgumentError normalize_model(
        ModelDefinition(;
            components=default_components(family),
            processes=default_processes(family),
            parameters=(invalid=1,),
        ),
    )
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
            shared_mortality=Parameter(ConstantDefault(0.1); axes=:plankton),
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

    @test_throws ArgumentError normalize_model(ModelDefinition(;
        components=(P=components.P, Z=components.Z),
        processes=(mortality=shared,),
        parameters=(shared_mortality=Parameter(
            ConstantDefault(0.1); axes=(:consumer, :prey)
        ),),
    ))

    accidental = ModelDefinition(;
        components=(P=components.P, Z=components.Z),
        processes=(
            mortality=Mortality(
                Agate.Processes.LinearMortality(); populations=(:P, :Z)
            ),
        ),
        parameters=(rate=Parameter(ConstantDefault(0.1); axes=:plankton),),
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
        parameters=(shared_mortality=Parameter(ConstantDefault(0.1); axes=:plankton),),
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
        parameters=(K=Parameter(ConstantDefault(0.1); axes=:plankton),),
    )
    @test_throws ArgumentError normalize_model(unknown_zero_slot)

    dependency_only = normalize_model(ModelDefinition(;
        components=(P=components.P,),
        processes=(mortality=Mortality(
            Agate.Processes.LinearMortality();
            populations=:P,
            bindings=(rate=:mortality_rate,),
        ),),
        parameters=(
            mortality_rate=Parameter(
                DerivedDefault(BindingDependencyDefault(); deps=(:middle,));
                axes=:plankton,
            ),
            middle=Parameter(DerivedDefault(BindingDependencyDefault(); deps=(:helper,))),
            helper=Parameter(ConstantDefault(1.0)),
        ),
    ))
    @test dependency_only.parameters.helper.spec.axes === nothing
    @test dependency_only.derived_parameter_order == (:middle, :mortality_rate)

    @test_throws ArgumentError normalize_model(ModelDefinition(;
        components=(P=components.P,),
        processes=(mortality=Mortality(
            Agate.Processes.LinearMortality();
            populations=:P,
            bindings=(rate=:mortality_rate,),
        ),),
        parameters=(
            mortality_rate=Parameter(ConstantDefault(0.1); axes=:plankton),
            unused=Parameter(ConstantDefault(1.0)),
        ),
    ))
end

@testset "Literal interaction defaults need no derivation dependencies" begin
    components = (
        P=Population(:carbon; size_structure=[1.0]),
        Z=Population(:carbon; size_structure=[10.0]),
    )
    processes = (
        grazing=Consumption(PreferentialGrazing(); consumers=:Z, resources=:P),
    )
    parameters = (
        maximum_rate=Parameter(ConstantDefault(1.0); axes=:plankton),
        half_saturation=Parameter(ConstantDefault(0.5); axes=:plankton),
        palatability=Parameter(ConstantDefault(1.0); axes=(:consumer, :prey)),
        assimilation=Parameter(ConstantDefault(0.7); axes=(:consumer, :prey)),
    )

    normalized = normalize_model(ModelDefinition(; components, processes, parameters))
    @test keys(normalized.parameters) ==
        (:maximum_rate, :half_saturation, :palatability, :assimilation)
    @test Tuple(def.spec.axes for def in normalized.parameters) ==
        (:plankton, :plankton, (:consumer, :prey), (:consumer, :prey))
end
