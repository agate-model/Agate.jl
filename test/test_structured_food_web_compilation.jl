using ForwardDiff

using Agate.Compilation:
    process_fluxes, group_fluxes, compile_tendencies, compile_model_tendencies
using Agate.Configuration:
    Population, Pool, realize_model_layout, component_tracers
using Agate.Construction: construct, define_tracer_functions
using Agate.Parameters: Parameter, NoDefault
using Agate.Processes:
    ModelDefinition, Growth, Light, NutrientResponse, Temperature, Consumption, Smith, Monod, Q10, HeterotrophicConsumption, PreferentialGrazing,
    normalize_model, participants, driver_identities, build_parameter_plan

function food_web_parameters()
    no_default(; axes=nothing) = Parameter(NoDefault(); axes)

    return (
        maximum_growth_rate=no_default(; axes=:plankton),
        alpha=no_default(; axes=:plankton),
        nutrient_half_saturation=no_default(; axes=:plankton),
        temperature_q10=no_default(),
        reference_temperature=no_default(),
        maximum_consumption_rate=no_default(; axes=:plankton),
        pom_half_saturation=no_default(),
        bacterial_assimilation=no_default(),
        maximum_predation_rate=no_default(; axes=:plankton),
        holling_half_saturation=no_default(; axes=:plankton),
        living_palatability_matrix=no_default(; axes=(:consumer, :prey)),
        living_assimilation_matrix=no_default(; axes=(:consumer, :prey)),
    )
end

function food_web_compilation(::Type{T}=Float64) where {T<:Real}
    components = (
        N=Pool(:nitrogen),
        D=Pool(:nitrogen),
        POM=Pool(:nitrogen; size_structure=T[0.5, 5]),
        P=Population(:nitrogen; size_structure=T[1]),
        B=Population(:nitrogen; size_structure=T[0.8]),
        M=Population(:nitrogen; size_structure=T[2]),
        Z=Population(:nitrogen; size_structure=T[10]),
    )
    temperature = Temperature(
        Q10(); bindings=(q10=:temperature_q10, reference_temperature=:reference_temperature)
    )
    processes = (
        growth_autotrophs=Growth(;
            populations=(:P, :M),
            factors=(
                temperature=temperature,
                nutrients=NutrientResponse(
                    Monod(); resource=:N, bindings=(K=:nutrient_half_saturation,)
                ),
                light=Light(
                    Smith(); driver=:PAR, bindings=(maximum_rate=:maximum_growth_rate,)
                ),
            ),
        ),
        consume_POM=Consumption(
            HeterotrophicConsumption();
            consumers=:B,
            resources=:POM,
            bindings=(
                maximum_rate=:maximum_consumption_rate,
                half_saturation=:pom_half_saturation,
                assimilation=:bacterial_assimilation,
            ),
            factors=(temperature=temperature,),
            unassimilated_products=:D,
        ),
        grazing_living=Consumption(
            PreferentialGrazing();
            consumers=(:M, :Z),
            resources=(:P, :B),
            bindings=(
                maximum_rate=:maximum_predation_rate,
                half_saturation=:holling_half_saturation,
                palatability=:living_palatability_matrix,
                assimilation=:living_assimilation_matrix,
            ),
            unassimilated_products=:D,
        ),
    )
    normalized = normalize_model(ModelDefinition(;
        components, processes, parameters=food_web_parameters()
    ))
    drivers = driver_identities(normalized)
    layout = realize_model_layout(
        components;
        scalar_type=T,
        interaction_roles=(consumers=(:M, :Z), prey=(:P, :B)),
        auxiliary_fields=drivers,
    )
    target_order = layout.tracer_order
    plan = build_parameter_plan(normalized, layout)
    compiled = compile_model_tendencies(normalized, layout, plan; target_order)
    return (; normalized, layout, plan, compiled, target_order)
end

function food_web_bgc(compilation)
    T = compilation.layout.scalar_type
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
        living_palatability_matrix=T[0.6 0.8; 0.7 0.9],
        living_assimilation_matrix=T[0.4 0.5; 0.35 0.45],
    )
    drivers = driver_identities(compilation.normalized)
    tracer_index = Agate.Runtime.build_tracer_index(compilation.layout)
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
    @test driver_identities(normalized) == (:PAR, :temperature)

    half_saturation = compilation.plan.parameters.pom_half_saturation
    assimilation = compilation.plan.parameters.bacterial_assimilation
    @test (half_saturation.storage_shape, half_saturation.storage_labels) ==
        ((2,), ((:POM_1, :POM_2),))
    @test (assimilation.storage_shape, assimilation.storage_labels) ==
        ((1, 2), ((:B_1,), (:POM_1, :POM_2)))

    consumption = normalized.processes.consume_POM
    fluxes = process_fluxes(
        consumption, normalized, layout, compilation.plan
    )
    growth_fluxes = process_fluxes(
        normalized.processes.growth_autotrophs, normalized, layout, compilation.plan
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
        X=Population(:nitrogen; size_structure=[0.4]),
        B=Population(:nitrogen; size_structure=[0.8]),
    )
    processes = (
        consume_POM=Consumption(
            HeterotrophicConsumption();
            consumers=:B,
            resources=:POM,
            bindings=(
                maximum_rate=:maximum_consumption_rate,
                half_saturation=:pom_half_saturation,
                assimilation=:bacterial_assimilation,
            ),
            unassimilated_products=:N,
        ),
    )
    parameters = (
        maximum_consumption_rate=Parameter(
            NoDefault();
            axes=:plankton,
        ),
        pom_half_saturation=Parameter(NoDefault()),
        bacterial_assimilation=Parameter(NoDefault()),
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
