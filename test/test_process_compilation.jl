using Agate.Compilation:
    TracerOp,
    ScalarParamOp,
    VecParamOp,
    MatParamOp,
    ComplementOp,
    process_fluxes,
    weight_sign
using Agate.Configuration: realize_components
using Agate.ModelFamilies: default_components
using Agate.Processes:
    ModelDefinition, MultiplicativeFactors, Smith, Monod, normalize_model

@testset "Process compilation" begin
    family = Agate.Models.NiPiZD.NiPiZDFamily()
    normalized = normalize_model(ModelDefinition(family))
    layout = realize_components(default_components(family))
    context = Agate.Configuration.parse_community(
        Float64,
        default_nipizd_community();
        biogeochem_tracers=(:N, :D),
        interaction_roles=(consumers=(:Z,), prey=(:P,)),
    )

    @testset "Growth" begin
        process = normalized.processes.growth_P
        fluxes = process_fluxes(process, normalized, layout, context)
        @test Tuple((flux.target, weight_sign(flux.weight)) for flux in fluxes) ==
            ((:P_1, 1), (:N, -1), (:P_2, 1), (:N, -1))

        rate = fluxes[1].rate
        @test rate.formulation isa MultiplicativeFactors
        @test rate.operands == (
            TracerOp{:P_1}(), VecParamOp{:maximum_growth_rate,3}()
        )
        @test Tuple(typeof(factor.formulation) for factor in rate.factors) == (Smith, Monod)
    end

    @testset "Living-resource consumption" begin
        process = normalized.processes.grazing_Z_on_P
        fluxes = process_fluxes(process, normalized, layout, context)
        @test Tuple((flux.target, weight_sign(flux.weight)) for flux in fluxes[1:3]) ==
            ((:P_1, -1), (:Z_1, 1), (:D, 1))
        @test fluxes[1].rate.operands == (
            TracerOp{:P_1}(),
            TracerOp{:Z_1}(),
            VecParamOp{:maximum_predation_rate,1}(),
            VecParamOp{:holling_half_saturation,1}(),
            MatParamOp{:palatability_matrix,1,1}(),
        )
        assimilation = MatParamOp{:assimilation_matrix,1,1}()
        @test fluxes[2].weight.operands == (assimilation,)
        @test fluxes[3].weight.operands == (ComplementOp(assimilation),)
    end

    @testset "Mortality" begin
        fluxes = ()
        for id in (
            :linear_mortality_P_to_N,
            :linear_mortality_P_to_D,
            :linear_mortality_Z_to_N,
            :quadratic_mortality_Z_to_D,
        )
            process = getproperty(normalized.processes, id)
            fluxes = (fluxes..., process_fluxes(process, normalized, layout, context)...)
        end
        @test Tuple(
            flux.rate.operands[1] for flux in fluxes if
            flux.process === :linear_mortality_P_to_N && weight_sign(flux.weight) == -1
        ) == (TracerOp{:P_1}(), TracerOp{:P_2}())
        @test fluxes[2].target === :N
        @test isempty(fluxes[2].weight.operands)
    end

    @testset "Remineralization" begin
        process = normalized.processes.remineralization_D
        fluxes = process_fluxes(process, normalized, layout, context)
        @test Tuple((flux.target, weight_sign(flux.weight)) for flux in fluxes) ==
            ((:D, -1), (:N, 1))
        @test fluxes[1].rate.operands ==
            (TracerOp{:D}(), ScalarParamOp{:detritus_remineralization}())
    end
end
