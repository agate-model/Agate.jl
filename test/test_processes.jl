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
    DirectRouting,
    ModelDefinition,
    ProductRouting,
    parameter_bindings,
    resolve_parameter_applicability,
    drivers,
    driver_identities,
    formulation,
    formulation_tag,
    factor_value,
    normalize_model,
    participants,
    process_rate,
    PreferentialGrazing

struct ExternalTestFormulation <: AbstractFormulation end

struct BindingDependencyDefault end
struct ExternalTestFactor{F<:AbstractFormulation} <: AbstractFactor
    formulation::F
end
Agate.Processes.formulation_tag(::ExternalTestFormulation) = :external_test

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
    @test formulation_tag(formulation(frank_nutrients)) === :frank_tnorm
    @test participants(growth) == (population=(:P,),)
    @test drivers(growth) == (light=:PAR,)

    external_factor = ExternalTestFactor(ExternalTestFormulation())
    external_model = normalize_model(ModelDefinition(;
        components=(P=Population(:nitrogen),),
        processes=(growth=Growth(; populations=:P, factors=(external=external_factor,)),),
    ))
    @test formulation_tag(formulation(external_factor)) === :external_test
    @test isempty(parameter_bindings(external_model))

    grazing = Consumption(
        PreferentialGrazing();
        consumers=:Z,
        resources=(:P, :B),
        routing=ProductRouting(DirectRouting(); destination=:D),
    )
    @test participants(grazing) == (consumer=(:Z,), resource=(:P, :B))
    @test formulation(grazing.routing) isa DirectRouting
    @test grazing.routing.retained === :D

    shared_driver_model = normalize_model(
        ModelDefinition(;
            components=(
                P=Population(:nitrogen), Z=Population(:nitrogen)
            ),
            processes=(
                growth_Z=Growth(;
                    populations=:Z, factors=(light=Light(Smith(); driver=:PAR),)
                ),
                growth_P=Growth(;
                    populations=:P, factors=(light=Light(Smith(); driver=:PAR),)
                ),
            ),
        ),
    )
    @test driver_identities(shared_driver_model) == (:PAR,)

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
            factors=(nutrients=NutrientResponse(Monod(); resource=:N),),
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

    # Normalization remains the trusted boundary when field constructors bypass authoring checks.
    bypassed_factor = ModelDefinition(;
        components=(P=Population(:nitrogen), N=Pool(:nitrogen)),
        processes=(growth=Growth(;
            populations=:P,
            factors=(
                light=Light(Monod(), :PAR, NamedTuple()),
                nutrients=NutrientResponse(Monod(); resource=:N),
            ),
        ),),
    )
    @test_throws ArgumentError normalize_model(bypassed_factor)

    bypassed_mortality = ModelDefinition(;
        components=(P=Population(:nitrogen),),
        processes=(mortality=Mortality(Monod(), (:P,), nothing, NamedTuple()),),
    )
    @test_throws ArgumentError normalize_model(bypassed_mortality)
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
    @test driver_identities(normalized) == (:PAR,)

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

    linear_P_to_N_rate = only(filter(parameter_bindings(normalized)) do binding
        binding.process === :linear_mortality_P_to_N && binding.slot === :rate
    end)
    linear_P_to_D_rate = only(filter(parameter_bindings(normalized)) do binding
        binding.process === :linear_mortality_P_to_D && binding.slot === :rate
    end)
    @test linear_P_to_N_rate.parameter === :linear_mortality
    @test linear_P_to_D_rate.parameter === :linear_detrital_mortality

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
        destinations=:D,
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
        binding.qualifier.population for binding in mortality_bindings
    ) == Set((:P, :Z))
    remineralization_binding = only(filter(
        binding -> binding.process === :remineralization,
        normalized.parameter_bindings,
    ))
    @test (
        remineralization_binding.parameter,
        remineralization_binding.axes,
        remineralization_binding.qualifier,
    ) == (:remineralization_rate, (), (source=:D,))

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
            remineralization=Remineralization(
                Agate.Processes.LinearRemineralization();
                sources=(:D, :E),
                destinations=:D,
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
                DerivedDefault(BindingDependencyDefault(); deps=(:helper,));
                axes=:plankton,
            ),
            helper=Parameter(ConstantDefault(1.0)),
        ),
    ))
    @test dependency_only.parameters.helper.spec.axes === nothing

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
