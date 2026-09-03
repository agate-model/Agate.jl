using Test
using Agate.Components: Plankton, Pool
using Agate.ModelFamilies: default_components, default_processes
using Agate.Parameters: ConstantDefault, DerivedDefault, ConstructionParameter, Parameter
using Agate.Processes:
    AbstractFactor, AbstractFormulation, FactorizedGrowth, Smith, Geider, Monod,
    NormalizedDroop, QuotaRegulatedMonod, Liebig, FrankTNorm, Q10, Growth, Light,
    NutrientLimitation, Temperature, NutrientResponse, QuotaResponse, NutrientUptake,
    FixedStoichiometry, Consumption, Mortality, Products, ModelDefinition,
    driver_identities, formulation, HeterotrophicConsumption, LinearMortality,
    QuadraticMortality, LinearRemineralization, canonicalize_model, participants,
    PreferentialGrazing, parameter_slots

import Agate.Processes: factor_inputs

struct ExternalTestFormulation <: AbstractFormulation end

struct BindingDependencyDefault end
struct MultiDriverTestFactor{F<:AbstractFormulation} <: AbstractFactor
    formulation::F
    bindings::NamedTuple
end
MultiDriverTestFactor(formulation) = MultiDriverTestFactor(formulation, (scale=:environment_scale,))
Agate.Processes.authored_parameter_bindings(factor::MultiDriverTestFactor) = factor.bindings
Agate.Processes.parameter_slots(::ExternalTestFormulation) =
    (Agate.Processes.ParameterSlot(:scale; domain=:positive),)
factor_inputs(::MultiDriverTestFactor) = (
    Agate.Processes.FactorDriver(:wind),
    Agate.Processes.FactorDriver(:temperature),
    Agate.Processes.FactorComponent(:P),
)
Agate.Processes.factor_value(::ExternalTestFormulation, wind, temperature, biomass, scale) =
    scale * (wind + 2temperature + 3biomass)

@testset "Process authoring and canonicalization" begin
    @test Agate.ModelDefinition === ModelDefinition

    light = Light(Smith(); driver=:PAR)
    response = NutrientResponse(Monod(); resource=:N)
    growth = Growth(;
        plankton=:P,
        reference_resource=:N,
        factors=(light=light, nutrients=response),
    )
    @test formulation(light) isa Smith
    @test formulation(response) isa Monod

    @test participants(growth) == (plankton=(:P,), resource=(:N,))

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
                P=Plankton(; states=(nitrogen=:nitrogen,), reference_state=:nitrogen),
                Z=Plankton(; states=(nitrogen=:nitrogen,), reference_state=:nitrogen),
                N=Pool(:nitrogen),
            ),
            processes=(
                growth_Z=Growth(;
                    plankton=:Z,
                    reference_resource=:N,
                    factors=(
                        light=Light(Smith(); driver=:PAR),
                        nutrients=NutrientResponse(Monod(); resource=:N),
                    ),
                ),
                growth_P=Growth(;
                    plankton=:P,
                    reference_resource=:N,
                    factors=(
                        light=Light(Smith(); driver=:PAR),
                        nutrients=NutrientResponse(Monod(); resource=:N),
                    ),
                ),
            ),
        ),
    )
    @test driver_identities(shared_driver_model) == (:PAR,)

    multi_driver_definition = ModelDefinition(;
        components=(P=Plankton(; states=(nitrogen=:nitrogen,), reference_state=:nitrogen, size_structure=[1.0]), N=Pool(:nitrogen)),
        processes=(growth=Growth(;
            plankton=:P, reference_resource=:N,
            factors=(environment=MultiDriverTestFactor(ExternalTestFormulation()),),
        ),),
        parameters=(maximum_rate=Parameter(1.0), environment_scale=Parameter(0.5)),
    )
    multi_driver_model = canonicalize_model(multi_driver_definition)
    @test driver_identities(multi_driver_model) == (:temperature, :wind)
    multi_driver_bgc = Agate.Construction.construct(multi_driver_definition)
    @test multi_driver_bgc(Val(:P_1), 0, 0, 0, 0, 2.0, 1.0, 2.0, 3.0) == 5.0

    invalid_growth = ModelDefinition(;
        components=(
            P=Plankton(; states=(nitrogen=:nitrogen,), reference_state=:nitrogen),
            N=Pool(:nitrogen),
        ),
        processes=(growth=Growth(;
            plankton=:P,
            reference_resource=:N,
            factors=(nutrients=NutrientResponse(Monod(); resource=:missing),),
        ),),
    )
    @test_throws ArgumentError canonicalize_model(invalid_growth)

    @test_throws MethodError Light(:smith; driver=:PAR)
    @test_throws ArgumentError Growth(;
        plankton=(), reference_resource=:N, factors=(light=Light(Smith(); driver=:PAR),)
    )
    @test_throws ArgumentError Consumption(PreferentialGrazing(); consumers=(), resources=:P)
    @test_throws ArgumentError Consumption(PreferentialGrazing(); consumers=:Z, resources=(:P, 1))

    for build_process in (
        () -> Growth(;
            plankton=(:P, :P),
            reference_resource=:N,
            factors=(light=Light(Smith(); driver=:PAR),)
        ),
        () -> Consumption(
            PreferentialGrazing(); consumers=(:Z, :Z), resources=:P
        ),
        () -> Mortality(Agate.Processes.LinearMortality(); plankton=(:P, :P)),
        () -> Agate.Processes.Remineralization(
            Agate.Processes.LinearRemineralization(); sources=(:D, :D), destination=:N
        ),
    )
        @test_throws ArgumentError build_process()
    end

    @test_nowarn canonicalize_model(ModelDefinition(;
        components=(
            P=Plankton(; states=(nitrogen=:nitrogen,), reference_state=:nitrogen),
            N=Pool(:nitrogen),
        ),
        processes=(growth=Growth(;
            plankton=:P,
            reference_resource=:N,
            factors=(light=Light(Smith(); driver=:PAR),),
        ),),
    ))

    wrong_element = ModelDefinition(;
        components=(
            P=Plankton(; states=(carbon=:carbon,), reference_state=:carbon),
            DIC=Pool(:carbon),
            DIN=Pool(:phosphorus),
        ),
        processes=(growth=Growth(;
            plankton=:P,
            reference_resource=:DIC,
            additional_resources=(nitrogen=:DIN,),
            factors=(
                light=Light(Geider(); driver=:PAR),
                nutrients=NutrientLimitation(
                    Liebig();
                    responses=(nitrogen=NutrientResponse(Monod(); resource=:DIN),),
                ),
            ),
            stoichiometry=FixedStoichiometry(; reference_element=:carbon),
        ),),
    )
    @test_throws ArgumentError canonicalize_model(wrong_element)

    # Invalid built-in formulation combinations are rejected by their concrete objects or factor contract.
    @test_throws MethodError Light(Monod(), :PAR, NamedTuple())
    @test_throws MethodError Mortality(Monod(), (:P,), nothing, NamedTuple())
    @test_throws ArgumentError canonicalize_model(ModelDefinition(;
        components=(
            P=Plankton(; states=(nitrogen=:nitrogen,), reference_state=:nitrogen),
            Z=Plankton(; states=(nitrogen=:nitrogen,), reference_state=:nitrogen),
        ),
        processes=(grazing=Consumption(PreferentialGrazing(); consumers=:Z, resources=:P, factors=(light=Light(Smith(); driver=:PAR),)),),
    ))
end

@testset "Built-in parameter domains" begin
    nodes = (
        FactorizedGrowth(), Smith(), Geider(), Monod(), NormalizedDroop(),
        QuotaRegulatedMonod(), FrankTNorm(), Q10(), PreferentialGrazing(),
        HeterotrophicConsumption(), LinearMortality(), QuadraticMortality(),
        LinearRemineralization(), Products((a=:A, b=:B); fractions=(a=:fraction_a,)),
        FixedStoichiometry(; reference_element=:carbon),
    )
    expected(name) = name in (:minimum_quota, :maximum_quota, :hill, :sharpness, :q10) ? :positive :
        name === :reference_temperature ? :finite :
        name in (:assimilation, :fraction) ? :unit_interval : :nonnegative
    @test all(slot.domain === expected(slot.name) for node in nodes for slot in parameter_slots(node))
    @test_throws ArgumentError Agate.Processes.ParameterSlot(:x, (:consumer, :consumer))
end

@testset "Canonicalization owns authored structure" begin
    single = Plankton(; states=(nitrogen=:nitrogen,), reference_state=:nitrogen)
    light = Light(Smith(); driver=:PAR)
    one_process(id, process, components) = ModelDefinition(;
        components, processes=NamedTuple{(id,)}((process,))
    )

    @test_nowarn canonicalize_model(one_process(
        :growth,
        Growth(;
            plankton=:P,
            reference_resource=:DIC,
            additional_resources=(nitrogen=:N,),
            stoichiometry=FixedStoichiometry(; reference_element=:carbon),
            factors=(light=light,),
        ),
        (
            P=Plankton(; states=(carbon=:carbon,), reference_state=:carbon),
            DIC=Pool(:carbon),
            N=Pool(:nitrogen),
        ),
    ))

    cases = (
        (
            "Growth additional resources require stoichiometry",
            one_process(:growth, Growth(;
                plankton=:P,
                reference_resource=:DIC,
                additional_resources=(nitrogen=:N,),
                factors=(light=light,),
            ), (
                P=Plankton(; states=(carbon=:carbon,), reference_state=:carbon),
                DIC=Pool(:carbon),
                N=Pool(:nitrogen),
            )),
            ("process :growth", "additional_resources", "FixedStoichiometry", "together"),
        ),
        (
            "Growth stoichiometry requires additional resources",
            one_process(:growth, Growth(;
                plankton=:P,
                reference_resource=:DIC,
                stoichiometry=FixedStoichiometry(; reference_element=:carbon),
                factors=(light=light,),
            ), (
                P=Plankton(; states=(carbon=:carbon,), reference_state=:carbon),
                DIC=Pool(:carbon),
            )),
            ("process :growth", "FixedStoichiometry", "additional_resources", "together"),
        ),
        (
            "Light factor process compatibility",
            one_process(
                :consume,
                Consumption(
                    PreferentialGrazing();
                    consumers=:Z, resources=:P, factors=(light=light,),
                ),
                (Z=single, P=single),
            ),
            ("process :consume", "Light", "Consumption", "not applicable"),
        ),
        (
            "Mortality participant type",
            one_process(:mortality, Mortality(
                Agate.Processes.LinearMortality(); plankton=:D
            ), (D=Pool(:nitrogen),)),
            ("process :mortality", "must be a Plankton"),
        ),
        (
            "Remineralization element",
            one_process(:remineralization, Agate.Processes.Remineralization(
                Agate.Processes.LinearRemineralization(); sources=:D, destination=:N
            ), (D=Pool(:nitrogen), N=Pool(:phosphorus))),
            ("process :remineralization", "remineralization source :D", "element"),
        ),
    )

    for (label, definition, fragments) in cases
        @testset "$label" begin
            message = canonicalization_error_message(definition)
            @test all(fragment -> occursin(fragment, message), fragments)
        end
    end

    parameter_cases = (
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
        P=Plankton(; states=(nitrogen=:nitrogen,), reference_state=:nitrogen, size_structure=[1.0]),
        Z=Plankton(; states=(nitrogen=:nitrogen,), reference_state=:nitrogen, size_structure=[10.0]),
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
        plankton=(:P, :Z),
        bindings=(rate=:shared_mortality,),
    )
    @test_throws MethodError Parameter(ConstantDefault(0.1); axes=:plankton)
    @test_throws ArgumentError ConstructionParameter(ConstantDefault(0.1); axes=:consumer)

    construction_bound_error = canonicalization_error_message(ModelDefinition(;
        components=(P=components.P, Z=components.Z),
        processes=(mortality=shared,),
        parameters=(shared_mortality=ConstructionParameter(ConstantDefault(0.1); axes=:plankton),),
    ))
    @test occursin("construction-only and cannot bind to a process slot", construction_bound_error)

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
                Agate.Processes.LinearMortality(); plankton=(:P, :Z)
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
    @test_throws ArgumentError canonicalize_model(ModelDefinition(;
        components=(D=components.D, E=components.E, R=components.R),
        processes=(remineralization=Agate.Processes.Remineralization(
            Agate.Processes.LinearRemineralization(); sources=(:D, :E), destination=:R,
            bindings=(rate=(D=:D_rate, E=:E_rate, extra=:D_rate),),
        ),),
        parameters=(D_rate=Parameter(0.1), E_rate=Parameter(0.2)),
    ))

    unknown = ModelDefinition(;
        components=(P=components.P,),
        processes=(
            mortality=Mortality(
                Agate.Processes.LinearMortality();
                plankton=:P,
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
                plankton=:P,
                reference_resource=:D,
                factors=(
                    light=Light(Smith(); driver=:PAR),
                    nutrients=NutrientLimitation(
                        Liebig();
                        responses=(nitrogen=NutrientResponse(Monod(); resource=:D),),
                        bindings=(missing=:shared_mortality,),
                    ),
                ),
            ),
        ),
        parameters=(half_saturation=Parameter(ConstantDefault(0.1)),),
    )
    @test_throws ArgumentError canonicalize_model(unknown_zero_slot)

    single_mortality = Mortality(
        Agate.Processes.LinearMortality();
        plankton=:P,
        bindings=(rate=:x,),
    )
    lifecycle_cases = (
        (
            (x=Parameter(0.1), unused=Parameter(1.0)),
            "use ConstructionParameter for construction-only",
        ),
        (
            (x=Parameter(0.1), helper=ConstructionParameter(1.0)),
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
