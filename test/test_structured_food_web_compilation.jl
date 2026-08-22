using ForwardDiff

using Agate.Compilation:
    TracerOp, VecParamOp, MatParamOp,
    ScalarParamOp, ComplementOp, process_fluxes, model_fluxes, group_fluxes,
    compile_tendencies, compile_model_tendencies, weight_sign
using Agate.Configuration:
    Population, Pool, realize_components, component_tracers, parse_community
using Agate.Construction: construct, define_tracer_functions
using Agate.Parameters:
    ParameterProvision, ParameterDefinition, NoDefault
using Agate.Processes:
    ModelDefinition, Growth, Light, NutrientResponse, Temperature, Consumption, Grazing,
    Smith, Monod, Q10, normalize_model, resolve_parameter_applicability, participants, factors,
    driver_identities

function food_web_parameters()
    slot(process, name) = ParameterProvision(process, name)
    no_default(
        name, shape, provides; runtime_path=(name,),
        axes=shape === :vector ? :plankton : nothing,
    ) = ParameterDefinition(name, NoDefault(); shape, axes, runtime_path, provides)

    return (
        no_default(:maximum_growth_rate, :vector,
            slot(:growth_autotrophs, :maximum_rate)),
        no_default(:alpha, :vector,
            slot(:growth_autotrophs, :alpha)),
        no_default(:nutrient_half_saturation, :vector,
            slot(:growth_autotrophs, :K)),
        no_default(:temperature_q10, :scalar, (
            slot(:growth_autotrophs, :q10),
            slot(:consume_POM, :q10),
        )),
        no_default(:reference_temperature, :scalar, (
            slot(:growth_autotrophs, :reference_temperature),
            slot(:consume_POM, :reference_temperature),
        )),
        no_default(:maximum_consumption_rate, :vector,
            slot(:consume_POM, :maximum_rate)),
        no_default(
            :pom_half_saturation, :vector, slot(:consume_POM, :half_saturation); axes=nothing
        ),
        no_default(:bacterial_assimilation, :matrix,
            slot(:consume_POM, :assimilation)),
        no_default(:maximum_predation_rate, :vector,
            slot(:grazing_living, :maximum_rate)),
        no_default(:holling_half_saturation, :vector,
            slot(:grazing_living, :half_saturation)),
        no_default(
            :living_palatability_matrix, :matrix,
            slot(:grazing_living, :palatability);
            runtime_path=(:interactions, :living_palatability),
            axes=(:consumer, :prey),
        ),
        no_default(
            :living_assimilation_matrix, :matrix,
            slot(:grazing_living, :assimilation);
            runtime_path=(:interactions, :living_assimilation),
            axes=(:consumer, :prey),
        ),
        no_default(:optimum_predator_prey_ratio, :vector,
            slot(:grazing_living, :optimum_predator_prey_ratio)),
        no_default(:specificity, :vector,
            slot(:grazing_living, :specificity)),
        no_default(:protection, :vector,
            slot(:grazing_living, :protection)),
        no_default(:assimilation_efficiency, :vector,
            slot(:grazing_living, :assimilation_efficiency)),
    )
end

function food_web_compilation(::Type{T}=Float64) where {T<:Real}
    components = (
        N=Pool(:nitrogen),
        D=Pool(:nitrogen),
        POM=Pool(:nitrogen; size_structure=T[0.5, 5]),
        P=Population(; currency=:nitrogen, size_structure=T[1]),
        B=Population(; currency=:nitrogen, size_structure=T[0.8]),
        M=Population(; currency=:nitrogen, size_structure=T[2]),
        Z=Population(; currency=:nitrogen, size_structure=T[10]),
    )
    temperature = Temperature(:q10)
    processes = (
        growth_autotrophs=Growth(;
            populations=(:P, :M),
            factors=(
                temperature=temperature,
                nutrients=NutrientResponse(:monod; resource=:N),
                light=Light(:smith; driver=:PAR),
            ),
        ),
        consume_POM=Consumption(
            :heterotrophic;
            consumer=:B,
            resource=:POM,
            unassimilated_destination=:D,
            factors=(temperature=temperature,),
        ),
        grazing_living=Grazing(
            :preferential;
            consumers=(:M, :Z),
            resources=(:P, :B),
            unassimilated_destination=:D,
        ),
    )
    normalized = normalize_model(ModelDefinition(;
        components, processes, parameters=food_web_parameters()
    ))
    community = (
        P=(n=1, diameters=T[1], pft=NamedTuple()),
        B=(n=1, diameters=T[0.8], pft=NamedTuple()),
        M=(n=1, diameters=T[2], pft=NamedTuple()),
        Z=(n=1, diameters=T[10], pft=NamedTuple()),
    )
    context = parse_community(
        T,
        community;
        biogeochem_tracers=(:N, :D, :POM_1, :POM_2),
        interaction_roles=(consumers=(:M, :Z), prey=(:P, :B)),
    )
    layout = realize_components(components; scalar_type=T)
    target_order = layout.tracer_order
    fluxes = model_fluxes(normalized, layout, context)
    compiled = compile_model_tendencies(normalized, layout, context; target_order)
    return (; normalized, layout, context, fluxes, compiled, target_order)
end

function food_web_bgc(compilation)
    T = compilation.context.scalar_type
    parameters = (
        maximum_growth_rate=T[2e-5, 0, 1.4e-5, 0],
        alpha=T[2e-6, 0, 1.6e-6, 0],
        nutrient_half_saturation=T[0.2, 0, 0.3, 0],
        temperature_q10=T(2),
        reference_temperature=T(20),
        maximum_consumption_rate=T[0, 1.5e-5, 0, 0],
        pom_half_saturation=T[0.15, 0.4],
        bacterial_assimilation=reshape(T[0.65, 0.75], 1, 2),
        maximum_predation_rate=T[0, 0, 6e-5, 9e-5],
        holling_half_saturation=T[1, 1, 0.12, 0.18],
        interactions=(
            living_palatability=T[0.6 0.8; 0.7 0.9],
            living_assimilation=T[0.4 0.5; 0.35 0.45],
        ),
    )
    drivers = driver_identities(compilation.normalized)
    tracer_index = Agate.Runtime.build_tracer_index(
        compilation.context,
        compilation.target_order,
        drivers;
        n_biogeochem_tracers=4,
    )
    factory = define_tracer_functions(
        parameters, compilation.compiled; auxiliary_fields=drivers, tracer_index
    )
    return factory(parameters)
end

function food_web_args(::Type{T}, temperature) where {T}
    return (
        zero(T), zero(T), zero(T), zero(T),
        T(5), T(0.1), T(0.5), T(0.2), T(0.05), T(0.03), T(0.02), T(0.04),
        T(100), T(temperature),
    )
end

@testset "Structured POM, bacteria, mixotrophy, and reusable factors" begin
    compilation = food_web_compilation()
    normalized = compilation.normalized
    layout = compilation.layout

    @test component_tracers(layout, :POM) == (:POM_1, :POM_2)
    @test participants(normalized.processes.consume_POM) == (
        consumer=(:B,), resource=(:POM,)
    )
    @test participants(normalized.processes.grazing_living).resource == (:P, :B)
    @test :POM ∉ participants(normalized.processes.grazing_living).resource
    @test :M ∈ participants(normalized.processes.growth_autotrophs).population
    @test :M ∈ participants(normalized.processes.grazing_living).consumer
    @test keys(factors(normalized.processes.growth_autotrophs)) ==
          (:light, :nutrients, :temperature)
    @test keys(factors(normalized.processes.consume_POM)) == (:temperature,)
    @test driver_identities(normalized) == (:PAR, :temperature)
    population_groups = (P=(:P,), B=(:B,), M=(:M,), Z=(:Z,))
    @test Agate.Construction._process_interaction_roles(normalized, population_groups) ==
          (consumers=(:M, :Z), prey=(:P, :B))

    applicability = resolve_parameter_applicability(normalized, layout)
    pom_K = only(a for a in applicability if a.binding.parameter === :pom_half_saturation)
    bacterial_assimilation = only(
        a for a in applicability if a.binding.parameter === :bacterial_assimilation
    )
    @test pom_K.axis_classes == ((:POM_1, :POM_2),)
    @test bacterial_assimilation.axis_classes == ((:B_1,), (:POM_1, :POM_2))

    consumption = normalized.processes.consume_POM
    fluxes = process_fluxes(
        consumption, normalized, layout, compilation.context
    )
    @test length(fluxes) == 6
    @test Tuple((flux.target, weight_sign(flux.weight)) for flux in fluxes) == (
        (:POM_1, -1), (:B_1, 1), (:D, 1),
        (:POM_2, -1), (:B_1, 1), (:D, 1),
    )
    @test fluxes[1].rate.operands == (
        TracerOp{:POM_1}(),
        TracerOp{:B_1}(),
        VecParamOp{:maximum_consumption_rate,2}(),
        VecParamOp{:pom_half_saturation,1}(),
    )
    assimilation = MatParamOp{:bacterial_assimilation,1,1}()
    @test fluxes[2].weight.operands == (assimilation,)
    @test fluxes[3].weight.operands == (ComplementOp(assimilation),)
    consumption_factor = only(fluxes[1].rate.factors)
    @test consumption_factor.formulation isa Q10
    @test consumption_factor.operands == (
        TracerOp{:temperature}(),
        ScalarParamOp{:temperature_q10}(),
        ScalarParamOp{:reference_temperature}(),
    )

    growth_fluxes = process_fluxes(
        normalized.processes.growth_autotrophs, normalized, layout, compilation.context
    )
    growth_rate = first(growth_fluxes).rate
    @test Tuple(typeof(factor.formulation) for factor in growth_rate.factors) ==
          (Smith, Monod, Q10)
    @test growth_rate.factors[3].operands == consumption_factor.operands

    @test length(compilation.fluxes) == 22
    grouped = group_fluxes(compilation.fluxes; target_order=compilation.target_order)
    @test map(length, grouped) == (
        N=2, D=6, POM_1=1, POM_2=1, P_1=3, B_1=4, M_1=3, Z_1=2
    )
    @test all(equation -> isbitstype(typeof(equation.f)), values(compilation.compiled))
    @test all(
        equation -> all(term -> isbitstype(typeof(term)), equation.f.terms),
        values(compilation.compiled),
    )

    bgc = food_web_bgc(compilation)
    args = food_web_args(Float64, 25)
    tendencies = map(target -> bgc(Val(target), args...), compilation.target_order)
    @test isapprox(sum(tendencies), 0; atol=10 * eps(sum(abs, tendencies)))

    consumption_grouped = group_fluxes(fluxes)
    consumption_compiled = compile_tendencies(consumption_grouped)
    growth_compiled = compile_tendencies(group_fluxes(growth_fluxes))
    args20 = food_web_args(Float64, 20)
    args30 = food_web_args(Float64, 30)
    @test process_compiler_isapprox(
        consumption_compiled.POM_1(bgc, args30...),
        2 * consumption_compiled.POM_1(bgc, args20...),
    )
    @test process_compiler_isapprox(
        growth_compiled.P_1(bgc, args30...),
        2 * growth_compiled.P_1(bgc, args20...),
    )

    derivative = ForwardDiff.derivative(0.5) do pom
        dynamic_args = (
            args[1:6]..., pom, args[8:end]...
        )
        consumption_compiled.POM_1(bgc, dynamic_args...)
    end
    @test isfinite(derivative)
    @test derivative < 0
end


@testset "Constructed consumer-resource storage axes" begin
    components = (
        N=Pool(:nitrogen),
        POM=Pool(:nitrogen; size_structure=[0.5, 1.0, 2.0]),
        X=Population(; currency=:nitrogen, size_structure=[0.4]),
        B=Population(; currency=:nitrogen, size_structure=[0.8]),
    )
    processes = (
        consume_POM=Consumption(
            :heterotrophic;
            consumer=:B,
            resource=:POM,
            unassimilated_destination=:N,
        ),
    )
    parameters = (
        ParameterDefinition(
            :maximum_consumption_rate,
            NoDefault();
            shape=:vector,
            axes=:plankton,
            provides=ParameterProvision(:consume_POM, :maximum_rate),
        ),
        ParameterDefinition(
            :pom_half_saturation,
            NoDefault();
            shape=:vector,
            provides=ParameterProvision(:consume_POM, :half_saturation),
        ),
        ParameterDefinition(
            :bacterial_assimilation,
            NoDefault();
            shape=:matrix,
            provides=ParameterProvision(:consume_POM, :assimilation),
        ),
    )
    definition = ModelDefinition(; components, processes, parameters)
    bgc = construct(
        definition;
        parameter_overrides=(
            maximum_consumption_rate=[100.0, 2.0],
            pom_half_saturation=[1.0, 3.0, 7.0],
            bacterial_assimilation=reshape([0.2, 0.4, 0.8], 1, 3),
        ),
    )

    @test bgc.parameters.maximum_consumption_rate == [100.0, 2.0]
    @test bgc.parameters.pom_half_saturation == [1.0, 3.0, 7.0]
    @test bgc.parameters.bacterial_assimilation == reshape([0.2, 0.4, 0.8], 1, 3)

    names = Agate.Introspection.tracer_names(bgc)
    state = (N=0.0, POM_1=1.0, POM_2=1.0, POM_3=1.0, X_1=5.0, B_1=1.0)
    args = (0.0, 0.0, 0.0, 0.0, Tuple(getproperty(state, name) for name in names)...)
    expected = (POM_1=-1.0, POM_2=-0.5, POM_3=-0.25, B_1=0.6, N=1.15, X_1=0.0)

    for (name, value) in pairs(expected)
        @test bgc(Val(name), args...) ≈ value
    end
    @test isapprox(sum(bgc(Val(name), args...) for name in names), 0.0; atol=1e-14)
end
