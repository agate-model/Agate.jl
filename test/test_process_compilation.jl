using Agate.Compilation: CompileContext, InputOp, ParameterOp, input_operand, process_fluxes
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

    context = CompileContext(normalized, layout, plan)
    growth_flux = first(process_fluxes(normalized.processes.growth_P, context))
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

struct ExtensionTransfer <: Agate.Processes.AbstractProcess
    source::Symbol
    destination::Symbol
end
struct ExtensionTransferRate <: Agate.Processes.AbstractFormulation end

Agate.Processes.formulation(::ExtensionTransfer) = ExtensionTransferRate()
Agate.Processes.participants(process::ExtensionTransfer) = (
    source=(process.source,), destination=(process.destination,)
)
Agate.Processes.normalize_process_facts(
    process::ExtensionTransfer, ::Symbol, ::NamedTuple
) = (; source=process.source, destination=process.destination)
Agate.Processes.process_rate(::ExtensionTransferRate, source) = source

function Agate.Compilation.process_fluxes(
    named::Agate.Processes.NamedProcess{P}, context::CompileContext
) where {P<:ExtensionTransfer}
    source = only(getproperty(context.layout.component_tracers, named.facts.source))
    destination = only(getproperty(context.layout.component_tracers, named.facts.destination))
    rate = Agate.Compilation.RateElement(
        ExtensionTransferRate(), (input_operand(context.layout, source),)
    )
    return (
        Agate.Compilation.FluxSpec(source, rate, Agate.Compilation.Weight{-1}()),
        Agate.Compilation.FluxSpec(destination, rate, Agate.Compilation.Weight{1}()),
    )
end

@testset "Internal custom process extension seam" begin
    definition = ModelDefinition(;
        components=(source=Agate.Configuration.Pool(:carbon), sink=Agate.Configuration.Pool(:carbon)),
        processes=(transfer=ExtensionTransfer(:source, :sink),),
    )
    bgc = Agate.Construction.construct(definition)
    @test bgc(Val(:source), 0, 0, 0, 0, 2.0, 1.0) == -2.0
    @test bgc(Val(:sink), 0, 0, 0, 0, 2.0, 1.0) == 2.0
end
