using Agate.Compilation:
    GrowthParameterBinding,
    GrowthBiomassContribution,
    GrowthResourceLossContribution,
    GrazingParameterBinding,
    GrazingResourceLossContribution,
    GrazingConsumerGainContribution,
    GrazingUnassimilatedContribution,
    MortalityParameterBinding,
    MortalityLossContribution,
    MortalityRoutingContribution,
    RemineralizationParameterBinding,
    RemineralizationSourceLossContribution,
    RemineralizationDestinationGainContribution,
    realize_process_topology,
    process_contributions,
    group_contributions
using Agate.Configuration: realize_components
using Agate.Factories: default_components
using Agate.Processes: ModelDefinition, normalize_model

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
        contributions = process_contributions(process, topology, binding)
        grouped = group_contributions(contributions; target_order=(:N, :P_1, :P_2))

        @test (binding.maximum_rate, binding.half_saturation, binding.alpha) ==
            (:maximum_growth_rate, :nutrient_half_saturation, :alpha)
        @test topology.population_tracers == (:P_1, :P_2)
        @test topology.population_indices == (3, 4)
        @test topology.resource_target === :N
        @test topology.light_driver === :PAR
        @test count(c -> c isa GrowthBiomassContribution, contributions) == 2
        @test count(c -> c isa GrowthResourceLossContribution, contributions) == 2
        @test map(length, grouped) == (N=2, P_1=1, P_2=1)
    end

    @testset "Grazing" begin
        process = normalized.processes.grazing_Z_on_P
        topology = realize_process_topology(process, layout, context)
        binding = GrazingParameterBinding(normalized, :grazing_Z_on_P)
        contributions = process_contributions(process, topology, binding)
        grouped = group_contributions(
            contributions; target_order=(:D, :Z_1, :Z_2, :P_1, :P_2)
        )

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
        @test count(c -> c isa GrazingResourceLossContribution, contributions) == 4
        @test count(c -> c isa GrazingConsumerGainContribution, contributions) == 4
        @test count(c -> c isa GrazingUnassimilatedContribution, contributions) == 4
        @test Tuple(
            (c.consumer_axis, c.resource_axis) for c in contributions if
            c isa GrazingResourceLossContribution
        ) == ((1, 1), (1, 2), (2, 1), (2, 2))
        @test map(length, grouped) == (D=4, Z_1=2, Z_2=2, P_1=2, P_2=2)
    end

    @testset "Mortality" begin
        contributions = ()
        for id in (:linear_mortality_P, :linear_mortality_Z, :quadratic_mortality_Z)
            process = getproperty(normalized.processes, id)
            topology = realize_process_topology(process, layout, context)
            binding = MortalityParameterBinding(normalized, id)
            contributions = (
                contributions..., process_contributions(process, topology, binding)...
            )
        end
        grouped = group_contributions(
            contributions; target_order=(:N, :D, :Z_1, :Z_2, :P_1, :P_2)
        )

        @test MortalityParameterBinding(normalized, :linear_mortality_P).rate ===
            :linear_mortality
        @test MortalityParameterBinding(normalized, :quadratic_mortality_Z).rate ===
            :quadratic_mortality
        @test MortalityParameterBinding(normalized, :linear_mortality_Z).routing_fraction ===
            :mortality_export_fraction
        @test count(c -> c isa MortalityLossContribution, contributions) == 6
        @test count(c -> c isa MortalityRoutingContribution, contributions) == 12
        @test Tuple(
            c.population_index for c in contributions if
            c isa MortalityLossContribution && c.process === :linear_mortality_P
        ) == (3, 4)
        @test map(length, grouped) == (N=6, D=6, Z_1=2, Z_2=2, P_1=1, P_2=1)
    end

    @testset "Remineralization" begin
        process = normalized.processes.remineralization_D
        topology = realize_process_topology(process, layout, context)
        binding = RemineralizationParameterBinding(normalized, :remineralization_D)
        contributions = process_contributions(process, topology, binding)

        @test binding.rate === :detritus_remineralization
        @test topology.source_tracers == (:D,)
        @test topology.destination_target === :N
        @test contributions[1] isa RemineralizationSourceLossContribution
        @test contributions[2] isa RemineralizationDestinationGainContribution
        @test (contributions[1].target, contributions[2].target) == (:D, :N)
    end
end
