using ForwardDiff

using Agate.Compilation: InputOp, ParameterOp, process_fluxes, weight_sign
using Agate.Construction: construct
using Agate.Library.Nutrients:
    monod_limitation, normalized_droop_limitation, quota_uptake_regulation
using Agate.Library.Photosynthesis: smith_light_limitation
using Agate.Processes: build_parameter_plan, driver_identities, normalize_model

quota_runtime_args(;
    DIC=10.0, DIN=1.0, PO4=1.0,
    P_1_carbon=1.0, P_1_nitrogen=0.1, P_1_phosphorus=0.01,
    P_2_carbon=0.0, P_2_nitrogen=0.0, P_2_phosphorus=0.0, PAR=100.0,
) = (0, 0, 0, 0, DIC, DIN, PO4, P_1_carbon, P_1_nitrogen, P_1_phosphorus,
     P_2_carbon, P_2_nitrogen, P_2_phosphorus, PAR)

function quota_tendencies(bgc, args)
    names = keys(bgc.equations)
    return NamedTuple{names}(Tuple(bgc(Val(name), args...) for name in names))
end

function quota_construct_error(overrides; definition=quota_definition())
    try
        construct(definition; parameter_overrides=overrides)
    catch error
        @test error isa ArgumentError
        return sprint(showerror, error)
    end
    @test false
    return ""
end

@testset "State-aware quota compilation" begin
    definition = quota_definition()
    normalized = normalize_model(definition)
    layout = Agate.Configuration.realize_components(
        definition.components; scalar_type=Float64,
        auxiliary_fields=driver_identities(normalized),
    )
    plan = build_parameter_plan(normalized, layout)
    growth = process_fluxes(normalized.processes.growth, normalized, layout, plan)
    uptake = process_fluxes(normalized.processes.nitrogen_uptake, normalized, layout, plan)

    @test Tuple((f.target, weight_sign(f.weight)) for f in growth) ==
        ((:P_1_carbon, 1), (:DIC, -1), (:P_2_carbon, 1), (:DIC, -1))
    @test Tuple((f.target, weight_sign(f.weight)) for f in uptake) ==
        ((:DIN, -1), (:P_1_nitrogen, 1), (:DIN, -1), (:P_2_nitrogen, 1))
    n_response = growth[1].rate.factors[2].operands[1].operands[1]
    @test n_response.operands == (
        InputOp{:P_1_nitrogen,5}(), InputOp{:P_1_carbon,4}(),
        ParameterOp{:minimum_nitrogen_quota,(1,)}(),
        ParameterOp{:maximum_nitrogen_quota,(1,)}(),
    )
    @test uptake[1].rate.operands[1:3] ==
        (InputOp{:DIN,2}(), InputOp{:P_1_nitrogen,5}(), InputOp{:P_1_carbon,4}())
    @test all(x -> isbitstype(typeof(x)), (growth[1].rate, uptake[1].rate))

    bgc = construct(definition)
    @test all(e -> isbitstype(typeof(e)), values(bgc.equations))
    @test all(e -> !hasproperty(e, :terms) || all(t -> isbitstype(typeof(t)), e.terms), values(bgc.equations))
    @test length(bgc.parameters.minimum_nitrogen_quota) == 2
    bgc32 = construct(definition; scalar_type=Float32)
    args32 = map(x -> x isa AbstractFloat ? Float32(x) : x, quota_runtime_args())
    @test bgc32(Val(:P_1_carbon), args32...) isa Float32

    tendency = quota_tendencies(bgc, quota_runtime_args())
    for (external, states) in (
        (:DIC, (:P_1_carbon, :P_2_carbon)),
        (:DIN, (:P_1_nitrogen, :P_2_nitrogen)),
        (:PO4, (:P_1_phosphorus, :P_2_phosphorus)),
    )
        transfer_values = (getproperty(tendency, external), (getproperty(tendency, s) for s in states)...)
        @test sum(transfer_values) == 0
        @test first(transfer_values) == -transfer_values[2]
    end
    @test all(>(0), (tendency.P_1_carbon, tendency.P_1_nitrogen, tendency.P_1_phosphorus))
    reference = 1.0
    for (state, internal, maximum_rate, K, minimum, maximum) in (
        (:P_1_nitrogen, 0.1, 0.1, 0.2, 0.05, 0.2),
        (:P_1_phosphorus, 0.01, 0.01, 0.02, 0.005, 0.02),
    )
        expected = maximum_rate * reference * monod_limitation(1.0, K) *
                   quota_uptake_regulation(internal, reference, minimum, maximum, 2.0)
        @test isapprox(getproperty(tendency, state), expected)
    end

    no_external = quota_tendencies(bgc, quota_runtime_args(; DIN=0.0, PO4=0.0))
    @test (no_external.P_1_nitrogen, no_external.P_1_phosphorus) == (0, 0)
    @test no_external.P_1_carbon > 0
    for (kwargs, state) in (
        ((P_1_nitrogen=0.05, P_1_phosphorus=0.015), :P_1_carbon),
        ((P_1_nitrogen=0.2,), :P_1_nitrogen),
    )
        @test getproperty(quota_tendencies(bgc, quota_runtime_args(; kwargs...)), state) == 0
    end

    expected_growth(nitrogen, phosphorus) = 0.5 * smith_light_limitation(100.0, 0.05, 0.5) *
        min(normalized_droop_limitation(nitrogen, 1.0, 0.05, 0.2),
            normalized_droop_limitation(phosphorus, 1.0, 0.005, 0.02))
    for (nitrogen, phosphorus) in ((0.075, 0.019), (0.19, 0.008))
        actual = quota_tendencies(
            bgc, quota_runtime_args(; P_1_nitrogen=nitrogen, P_1_phosphorus=phosphorus)
        ).P_1_carbon
        @test isapprox(actual, expected_growth(nitrogen, phosphorus))
    end

    zero = quota_tendencies(bgc, quota_runtime_args(;
        P_1_carbon=0.0, P_1_nitrogen=0.1, P_1_phosphorus=0.01,
    ))
    @test all(isfinite, values(zero)) && all(iszero, values(zero))
    carbon_growth(nitrogen) = bgc(Val(:P_1_carbon), quota_runtime_args(;
        P_1_nitrogen=nitrogen, P_1_phosphorus=0.018,
    )...)
    derivative = ForwardDiff.derivative(carbon_growth, 0.075)
    @test isfinite(derivative) && derivative > 0
end

@testset "Realized quota parameter validation" begin
    cases = (
        ((minimum_nitrogen_quota=[0.0, 0.05],), :growth, :minimum_nitrogen_quota, "0.0"),
        ((maximum_nitrogen_quota=[0.04, 0.2],), :growth, :maximum_nitrogen_quota, "0.04"),
        ((nitrogen_uptake_hill=[0.0, 2.0],), :nitrogen_uptake, :nitrogen_uptake_hill, "0.0"),
        ((nitrogen_half_saturation=[-0.1, 0.2],), :nitrogen_uptake, :nitrogen_half_saturation, "-0.1"),
        ((maximum_nitrogen_uptake=[-0.1, 0.1],), :nitrogen_uptake, :maximum_nitrogen_uptake, "-0.1"),
    )
    for (overrides, process, parameter, value) in cases
        message = quota_construct_error(overrides)
        @test all(occursin(fragment, message) for fragment in
            ("process :$process", "path", "P_1", String(parameter), value))
    end

    global_definition = Agate.Processes.ModelDefinition(;
        components=merge(quota_components(), (
            Q=Agate.Configuration.Population(:carbon; size_structure=[3.0]),
        )),
        processes=quota_processes(), parameters=quota_parameters(),
    )
    global_bgc = construct(global_definition; parameter_overrides=(
        nitrogen_uptake_hill=[2.0, 2.0, -1.0],
    ))
    @test global_bgc.parameters.nitrogen_uptake_hill == [2.0, 2.0, -1.0]
end
