using Agate.Compilation: CompileContext, input_operand
using Agate.Processes: ModelDefinition
using Agate.Configuration: component_tracers

struct ExtensionTransfer <: Agate.Processes.AbstractProcess
    source::Symbol
    destination::Symbol
end
struct ExtensionTransferRate <: Agate.Processes.AbstractFormulation end

Agate.Processes.formulation(::ExtensionTransfer) = ExtensionTransferRate()
Agate.Processes.participants(process::ExtensionTransfer) = (
    source=(process.source,), destination=(process.destination,)
)
Agate.Processes.process_facts(
    process::ExtensionTransfer, ::Symbol, ::NamedTuple
) = (; source=process.source, destination=process.destination)
Agate.Processes.process_rate(::ExtensionTransferRate, source) = source

function Agate.Compilation.process_fluxes(
    named::Agate.Processes.NamedProcess{P}, context::CompileContext
) where {P<:ExtensionTransfer}
    source = only(component_tracers(context.layout, named.facts.source))
    destination = only(component_tracers(context.layout, named.facts.destination))
    rate = Agate.Compilation.RateElement(
        ExtensionTransferRate(), (input_operand(context.layout, source),)
    )
    return (
        Agate.Compilation.FluxSpec(source, rate, Agate.Compilation.Weight{-1}()),
        Agate.Compilation.FluxSpec(destination, rate, Agate.Compilation.Weight{1}()),
    )
end

@testset "Custom process extension" begin
    definition = ModelDefinition(;
        components=(source=Agate.Configuration.Pool(:carbon), sink=Agate.Configuration.Pool(:carbon)),
        processes=(transfer=ExtensionTransfer(:source, :sink),),
    )
    bgc = Agate.Construction.construct(definition)
    @test bgc(Val(:source), 0, 0, 0, 0, 2.0, 1.0) == -2.0
    @test bgc(Val(:sink), 0, 0, 0, 0, 2.0, 1.0) == 2.0
end
