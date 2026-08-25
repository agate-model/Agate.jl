using Agate.Compilation:
    TracerOp,
    ScalarParamOp,
    VecParamOp,
    MatParamOp,
    ComplementOp,
    OneMinusSumOp,
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
        process = normalized.processes.linear_mortality_P
        fluxes = process_fluxes(process, normalized, layout, context)

        @test Tuple((flux.target, weight_sign(flux.weight)) for flux in fluxes) == (
            (:P_1, -1), (:D, 1), (:N, 1),
            (:P_2, -1), (:D, 1), (:N, 1),
        )
        @test Tuple(flux.rate.operands[1] for flux in fluxes if weight_sign(flux.weight) == -1) ==
            (TracerOp{:P_1}(), TracerOp{:P_2}())
        @test fluxes[2].weight.operands ==
            (ComplementOp(ScalarParamOp{:mortality_export_fraction}()),)
        @test fluxes[3].weight.operands ==
            (ScalarParamOp{:mortality_export_fraction}(),)
    end

    @testset "Arbitrary product partition" begin
        components = (
            P=default_components(family).P,
            A=Agate.Configuration.Pool(:nitrogen),
            B=Agate.Configuration.Pool(:nitrogen),
            R=Agate.Configuration.Pool(:nitrogen),
        )
        process = Agate.Processes.Mortality(
            Agate.Processes.LinearMortality();
            populations=:P,
            bindings=(rate=:linear_mortality,),
            products=Agate.Processes.Products(
                (a=:A, b=:B, c=:R);
                fractions=(a=:fraction_a, b=:fraction_b),
            ),
        )
        family_parameters = Agate.Parameters.parameter_definitions(family)
        parameters = (
            linear_mortality=family_parameters.linear_mortality,
            fraction_a=Agate.Parameters.Parameter(Agate.Parameters.ConstantDefault(0.2)),
            fraction_b=Agate.Parameters.Parameter(Agate.Parameters.ConstantDefault(0.3)),
        )
        normalized_partition = normalize_model(ModelDefinition(;
            components, processes=(partition=process,), parameters,
        ))
        fluxes = process_fluxes(
            normalized_partition.processes.partition,
            normalized_partition,
            realize_components(components),
            context,
        )
        products = fluxes[2:4]
        @test Tuple(flux.target for flux in products) == (:A, :B, :R)
        @test products[1].weight.operands == (ScalarParamOp{:fraction_a}(),)
        @test products[2].weight.operands == (ScalarParamOp{:fraction_b}(),)
        @test products[3].weight.operands == (
            OneMinusSumOp((ScalarParamOp{:fraction_a}(), ScalarParamOp{:fraction_b}())),
        )
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
