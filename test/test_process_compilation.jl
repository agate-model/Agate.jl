using Agate.Compilation: InputOp, ParameterOp, input_operand, process_fluxes
using Agate.Processes:
    ModelDefinition, driver_identities, normalize_model, build_parameter_plan

@testset "Static compiler IR contract" begin
    family = Agate.Models.NiPiZD.NiPiZDFamily()
    normalized = normalize_model(ModelDefinition(family))
    layout = default_nipizd_layout(; auxiliary_fields=driver_identities(normalized))
    plan = build_parameter_plan(normalized, layout)

    input = input_operand(layout, :P_1)
    @test input isa InputOp
    @test only(typeof(input).parameters) isa Int
    @test only(typeof(input_operand(layout, :PAR)).parameters) == length(layout.tracer_order) + 1

    growth_flux = first(process_fluxes(
        normalized.processes.growth_P, normalized, layout, plan
    ))
    parameter = only(op for op in growth_flux.rate.operands if op isa ParameterOp)
    parameter_indices = typeof(parameter).parameters[2]
    @test parameter_indices isa Tuple
    @test !isempty(parameter_indices)
    @test all(index -> index isa Int, parameter_indices)

    bgc = Agate.Models.NiPiZD.construct()
    @test all(equation -> isbitstype(typeof(equation)), values(bgc.equations))
    @test all(
        equation -> !hasproperty(equation, :terms) ||
                    all(term -> isbitstype(typeof(term)), equation.terms),
        values(bgc.equations),
    )
    @test !any(type -> type === Any, fieldtypes(typeof(bgc.equations)))

    @test_throws ArgumentError input_operand(layout, :temperature)
end
