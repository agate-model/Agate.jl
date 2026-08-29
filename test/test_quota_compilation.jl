using ForwardDiff

using Agate.Construction: construct
using Agate.Library.Nutrients:
    monod_limitation, normalized_droop_limitation, quota_uptake_regulation
using Agate.Library.Photosynthesis: smith_light_limitation

@testset "State-aware quota compilation" begin
    definition = quota_definition()
    bgc = construct(definition)
    args = quota_runtime_args()
    @test isfinite(@inferred(bgc(Val(:P_1_carbon), args...)))
    @test length(bgc.parameters.minimum_nitrogen_quota) == 2
    bgc32 = construct(definition; scalar_type=Float32)
    args32 = map(x -> x isa AbstractFloat ? Float32(x) : x, quota_runtime_args())
    @test bgc32(Val(:P_1_carbon), args32...) isa Float32

    tendency = model_tendencies(bgc, args)
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
    for (state, internal, maximum_rate, half_saturation, minimum, maximum) in (
        (:P_1_nitrogen, 0.1, 0.1, 0.2, 0.05, 0.2),
        (:P_1_phosphorus, 0.01, 0.01, 0.02, 0.005, 0.02),
    )
        expected = maximum_rate * reference * monod_limitation(1.0, half_saturation) *
                   quota_uptake_regulation(internal, reference, minimum, maximum, 2.0)
        @test isapprox(getproperty(tendency, state), expected)
    end

    no_external = model_tendencies(bgc, quota_runtime_args(; DIN=0.0, PO4=0.0))
    @test (no_external.P_1_nitrogen, no_external.P_1_phosphorus) == (0, 0)
    @test no_external.P_1_carbon > 0
    for (kwargs, state) in (
        ((P_1_nitrogen=0.05, P_1_phosphorus=0.015), :P_1_carbon),
        ((P_1_nitrogen=0.2,), :P_1_nitrogen),
    )
        @test getproperty(model_tendencies(bgc, quota_runtime_args(; kwargs...)), state) == 0
    end

    expected_growth(nitrogen, phosphorus) = 0.5 * smith_light_limitation(100.0, 0.05, 0.5) *
        min(normalized_droop_limitation(nitrogen, 1.0, 0.05, 0.2),
            normalized_droop_limitation(phosphorus, 1.0, 0.005, 0.02))
    for (nitrogen, phosphorus) in ((0.075, 0.019), (0.19, 0.008))
        actual = model_tendencies(
            bgc, quota_runtime_args(; P_1_nitrogen=nitrogen, P_1_phosphorus=phosphorus)
        ).P_1_carbon
        @test isapprox(actual, expected_growth(nitrogen, phosphorus))
    end

    zero_biomass = model_tendencies(bgc, quota_runtime_args(;
        P_1_carbon=0.0, P_1_nitrogen=0.1, P_1_phosphorus=0.01,
    ))
    @test all(isfinite, values(zero_biomass)) && all(iszero, values(zero_biomass))
    carbon_growth(nitrogen) = bgc(Val(:P_1_carbon), quota_runtime_args(;
        P_1_nitrogen=nitrogen, P_1_phosphorus=0.018,
    )...)
    derivative = ForwardDiff.derivative(carbon_growth, 0.075)
    @test isfinite(derivative) && derivative > 0
end


@testset "Hybrid explicit and implicit elemental growth" begin
    components = merge(quota_components(), (
        P=Agate.Configuration.Plankton(;
            states=(:carbon, :nitrogen), reference_state=:carbon
        ),
    ))
    quota = quota_processes()
    growth = Agate.Processes.Growth(;
        plankton=:P,
        reference_resource=:DIC,
        additional_resources=(phosphorus=:PO4,),
        stoichiometry=Agate.Processes.FixedStoichiometry(;
            reference_element=:carbon,
            bindings=(ratio=(phosphorus=:phosphorus_to_carbon,),),
        ),
        bindings=quota.growth.bindings,
        factors=(nutrients=Agate.Processes.NutrientLimitation(
            Agate.Processes.Liebig();
            responses=(
                nitrogen=quota.growth.factors.nutrients.responses.nitrogen,
                phosphorus=Agate.Processes.NutrientResponse(
                    Agate.Processes.Monod(); resource=:PO4,
                    bindings=(half_saturation=:phosphorus_half_saturation,),
                ),
            ),
        ),),
    )
    parameters = (
        maximum_growth_rate=Agate.Parameters.Parameter(0.5),
        phosphorus_to_carbon=Agate.Parameters.Parameter(0.1),
        phosphorus_half_saturation=Agate.Parameters.Parameter(0.0),
        minimum_nitrogen_quota=Agate.Parameters.Parameter(0.1),
        maximum_nitrogen_quota=Agate.Parameters.Parameter(0.2),
        maximum_nitrogen_uptake=Agate.Parameters.Parameter(0.2),
        nitrogen_half_saturation=Agate.Parameters.Parameter(0.0),
        nitrogen_uptake_hill=Agate.Parameters.Parameter(1.0),
    )
    bgc = construct(Agate.Processes.ModelDefinition(;
        components, processes=(growth=growth, nitrogen_uptake=quota.nitrogen_uptake), parameters
    ))

    growth_rate, uptake_rate = 2 / 3, 0.2
    test_tendencies(
        bgc,
        (DIC=10.0, DIN=10.0, PO4=10.0, P_carbon=2.0, P_nitrogen=0.3),
        (
            DIC=-growth_rate, DIN=-uptake_rate, PO4=-0.1 * growth_rate,
            P_carbon=growth_rate, P_nitrogen=uptake_rate,
        ),
    )
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
        message = argument_error_message(() -> construct(
            quota_definition(); parameter_overrides=overrides
        ))
        @test all(occursin(fragment, message) for fragment in
            ("process :$process", "path", "P_1", String(parameter), value))
    end

end
