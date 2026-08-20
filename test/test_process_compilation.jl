using Agate.Compilation:
    TracerOp,
    ClassOp,
    ScalarParamOp,
    VecParamOp,
    InteractionParamOp,
    ComplementOp,
    GrowthParameterBinding,
    GrazingParameterBinding,
    MortalityParameterBinding,
    RemineralizationParameterBinding,
    realize_process_topology,
    process_fluxes,
    group_fluxes,
    weight_sign
using Agate.Configuration: realize_components
using Agate.Factories: default_components
using Agate.Processes: ModelDefinition, Smith, Monod, normalize_model

@testset "Process compilation" begin
    factory = Agate.Models.NiPiZD.NiPiZDFactory()
    normalized = normalize_model(ModelDefinition(factory))
    layout = realize_components(default_components(factory))
    context = Agate.Configuration.parse_community(
        Float64, default_nipizd_community(); biogeochem_tracers=(:N, :D)
    )

    @testset "Growth" begin
        process = normalized.processes.growth_P
        topology = realize_process_topology(process, layout, context)
        binding = GrowthParameterBinding(normalized, :growth_P)
        fluxes = process_fluxes(process, topology, binding)
        grouped = group_fluxes(fluxes; target_order=(:N, :P_1, :P_2))

        @test (binding.maximum_rate, binding.half_saturation, binding.alpha) ==
            (:maximum_growth_rate, :nutrient_half_saturation, :alpha)
        @test topology.population_tracers == (:P_1, :P_2)
        @test topology.population_indices == (3, 4)
        @test topology.resource_target === :N
        @test topology.light_driver === :PAR
        @test Tuple((flux.target, weight_sign(flux.weight)) for flux in fluxes) ==
            ((:P_1, 1), (:N, -1), (:P_2, 1), (:N, -1))
        @test typeof.(fluxes[1].rate.formulation) == (Smith, Monod)
        @test fluxes[1].rate.operands == (
            ClassOp{3}(),
            TracerOp{:N}(),
            TracerOp{:PAR}(),
            VecParamOp{:maximum_growth_rate,3}(),
            VecParamOp{:nutrient_half_saturation,3}(),
            VecParamOp{:alpha,3}(),
        )
        @test map(length, grouped) == (N=2, P_1=1, P_2=1)
    end

    @testset "Grazing" begin
        process = normalized.processes.grazing_Z_on_P
        topology = realize_process_topology(process, layout, context)
        binding = GrazingParameterBinding(normalized, :grazing_Z_on_P)
        fluxes = process_fluxes(process, topology, binding)
        grouped = group_fluxes(fluxes; target_order=(:D, :Z_1, :Z_2, :P_1, :P_2))

        @test (
            binding.maximum_rate,
            binding.half_saturation,
            binding.palatability,
            binding.assimilation,
        ) == (
            :maximum_predation_rate,
            :holling_half_saturation,
            :palatability_matrix,
            :assimilation_matrix,
        )
        @test topology.consumer_tracers == (:Z_1, :Z_2)
        @test topology.consumer_indices == (1, 2)
        @test topology.resource_tracers == (:P_1, :P_2)
        @test topology.resource_indices == (3, 4)
        @test topology.unassimilated_target === :D
        @test length(fluxes) == 12
        @test Tuple((flux.target, weight_sign(flux.weight)) for flux in fluxes[1:3]) ==
            ((:P_1, -1), (:Z_1, 1), (:D, 1))
        @test fluxes[1].rate.operands == (
            ClassOp{3}(),
            ClassOp{1}(),
            VecParamOp{:maximum_predation_rate,1}(),
            VecParamOp{:holling_half_saturation,1}(),
            InteractionParamOp{:palatability,1,1}(),
        )
        assimilation = InteractionParamOp{:assimilation,1,1}()
        @test fluxes[2].weight.operands == (assimilation,)
        @test fluxes[3].weight.operands == (ComplementOp(assimilation),)
        @test Tuple(flux.rate.operands[5] for flux in fluxes[1:3:end]) == (
            InteractionParamOp{:palatability,1,1}(),
            InteractionParamOp{:palatability,1,2}(),
            InteractionParamOp{:palatability,2,1}(),
            InteractionParamOp{:palatability,2,2}(),
        )
        @test map(length, grouped) == (D=4, Z_1=2, Z_2=2, P_1=2, P_2=2)
    end

    @testset "Mortality" begin
        fluxes = ()
        for id in (:linear_mortality_P, :linear_mortality_Z, :quadratic_mortality_Z)
            process = getproperty(normalized.processes, id)
            topology = realize_process_topology(process, layout, context)
            binding = MortalityParameterBinding(normalized, id)
            fluxes = (fluxes..., process_fluxes(process, topology, binding)...)
        end
        grouped = group_fluxes(
            fluxes; target_order=(:N, :D, :Z_1, :Z_2, :P_1, :P_2)
        )

        @test MortalityParameterBinding(normalized, :linear_mortality_P).rate ===
            :linear_mortality
        @test MortalityParameterBinding(normalized, :quadratic_mortality_Z).rate ===
            :quadratic_mortality
        @test MortalityParameterBinding(normalized, :linear_mortality_Z).routing_fraction ===
            :mortality_export_fraction
        @test length(fluxes) == 18
        @test Tuple(
            flux.rate.operands[1] for flux in fluxes if
            flux.process === :linear_mortality_P && weight_sign(flux.weight) == -1
        ) == (ClassOp{3}(), ClassOp{4}())
        fraction = ScalarParamOp{:mortality_export_fraction}()
        @test fluxes[2].weight.operands == (ComplementOp(fraction),)
        @test fluxes[3].weight.operands == (fraction,)
        @test map(length, grouped) == (N=6, D=6, Z_1=2, Z_2=2, P_1=1, P_2=1)
    end

    @testset "Remineralization" begin
        process = normalized.processes.remineralization_D
        topology = realize_process_topology(process, layout, context)
        binding = RemineralizationParameterBinding(normalized, :remineralization_D)
        fluxes = process_fluxes(process, topology, binding)

        @test binding.rate === :detritus_remineralization
        @test topology.source_tracers == (:D,)
        @test topology.destination_target === :N
        @test Tuple((flux.target, weight_sign(flux.weight)) for flux in fluxes) ==
            ((:D, -1), (:N, 1))
        @test fluxes[1].rate.operands ==
            (TracerOp{:D}(), ScalarParamOp{:detritus_remineralization}())
    end
end
