using Agate.Compilation:
    InputOp,
    ParameterOp,
    ComplementOp,
    OneMinusSumOp,
    process_fluxes,
    weight_sign
using Agate.Configuration: realize_model_layout
using Agate.ModelFamilies: default_components
using Agate.Processes:
    ModelDefinition, MultiplicativeFactors, Smith, Monod, driver_identities, normalize_model, build_parameter_plan

@testset "Process compilation" begin
    family = Agate.Models.NiPiZD.NiPiZDFamily()
    normalized = normalize_model(ModelDefinition(family))
    layout = default_nipizd_layout(; auxiliary_fields=driver_identities(normalized))
    plan = build_parameter_plan(normalized, layout)

    @testset "Growth" begin
        process = normalized.processes.growth_P
        fluxes = process_fluxes(process, normalized, layout, plan)
        @test Tuple((flux.target, weight_sign(flux.weight)) for flux in fluxes) ==
            ((:P_1, 1), (:N, -1), (:P_2, 1), (:N, -1))

        rate = fluxes[1].rate
        @test rate.formulation isa MultiplicativeFactors
        @test rate.operands == (
            InputOp{:P_1,5}(), ParameterOp{:maximum_growth_rate,(3,)}()
        )
        @test Tuple(typeof(factor.formulation) for factor in rate.factors) == (Smith, Monod)
        @test rate.factors[1].operands[1] == InputOp{:PAR,7}()
    end

    @testset "Living-resource consumption" begin
        process = normalized.processes.grazing_Z_on_P
        fluxes = process_fluxes(process, normalized, layout, plan)
        @test Tuple((flux.target, weight_sign(flux.weight)) for flux in fluxes[1:3]) ==
            ((:P_1, -1), (:Z_1, 1), (:D, 1))
        @test fluxes[1].rate.operands == (
            InputOp{:P_1,5}(),
            InputOp{:Z_1,3}(),
            ParameterOp{:maximum_predation_rate,(1,)}(),
            ParameterOp{:holling_half_saturation,(1,)}(),
            ParameterOp{:palatability_matrix,(1,1)}(),
        )
        assimilation = ParameterOp{:assimilation_matrix,(1,1)}()
        @test fluxes[2].weight.operands == (assimilation,)
        @test fluxes[3].weight.operands == (ComplementOp(assimilation),)
    end

    @testset "Mortality" begin
        process = normalized.processes.linear_mortality_P
        fluxes = process_fluxes(process, normalized, layout, plan)

        @test Tuple((flux.target, weight_sign(flux.weight)) for flux in fluxes) == (
            (:P_1, -1), (:D, 1), (:N, 1),
            (:P_2, -1), (:D, 1), (:N, 1),
        )
        @test Tuple(flux.rate.operands[1] for flux in fluxes if weight_sign(flux.weight) == -1) ==
            (InputOp{:P_1,5}(), InputOp{:P_2,6}())
        @test fluxes[2].weight.operands ==
            (ComplementOp(ParameterOp{:mortality_export_fraction,()}()),)
        @test fluxes[3].weight.operands ==
            (ParameterOp{:mortality_export_fraction,()}(),)
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
        realization = default_nipizd_realization()
        partition_layout = realize_model_layout(
            components;
            population_groups=(P=(:P,),),
            group_diameters=(P=realization.group_diameters.P,),
        )
        partition_plan = build_parameter_plan(normalized_partition, partition_layout)
        fluxes = process_fluxes(
            normalized_partition.processes.partition, normalized_partition, partition_layout, partition_plan
        )
        products = fluxes[2:4]
        @test Tuple(flux.target for flux in products) == (:A, :B, :R)
        @test products[1].weight.operands == (ParameterOp{:fraction_a,()}(),)
        @test products[2].weight.operands == (ParameterOp{:fraction_b,()}(),)
        @test products[3].weight.operands == (
            OneMinusSumOp((ParameterOp{:fraction_a,()}(), ParameterOp{:fraction_b,()}())),
        )
    end

    @testset "Remineralization" begin
        process = normalized.processes.remineralization_D
        fluxes = process_fluxes(process, normalized, layout, plan)
        @test Tuple((flux.target, weight_sign(flux.weight)) for flux in fluxes) ==
            ((:D, -1), (:N, 1))
        @test fluxes[1].rate.operands ==
            (InputOp{:D,2}(), ParameterOp{:detritus_remineralization,()}())
    end
end
