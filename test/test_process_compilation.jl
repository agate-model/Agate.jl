using Agate.Compilation: CompileContext, input_operand, process_parameter_operands
using Agate.Processes: ModelDefinition
using Agate.Parameters: Parameter
using Agate.Configuration: component_tracers

struct ExtensionTransfer <: Agate.Processes.AbstractProcess
    source::Symbol
    destination::Symbol
    bindings::NamedTuple
end
struct ExtensionTransferRate <: Agate.Processes.AbstractFormulation end

ExtensionTransfer(source::Symbol, destination::Symbol; bindings::NamedTuple) =
    ExtensionTransfer(source, destination, bindings)

Agate.Processes.formulation(::ExtensionTransfer) = ExtensionTransferRate()
Agate.Processes.authored_parameter_bindings(process::ExtensionTransfer) = process.bindings
Agate.Processes.parameter_slots(::ExtensionTransferRate) =
    (Agate.Processes.ParameterSlot(:rate),)
Agate.Processes.participants(process::ExtensionTransfer) = (
    source=(process.source,), destination=(process.destination,)
)
Agate.Processes.process_facts(
    process::ExtensionTransfer, ::Symbol, ::NamedTuple
) = (; source=process.source, destination=process.destination)
Agate.Processes.process_rate(::ExtensionTransferRate, source, rate) = source * rate

function Agate.Compilation.process_fluxes(
    named::Agate.Processes.CanonicalProcess{P}, context::CompileContext
) where {P<:ExtensionTransfer}
    source = only(component_tracers(context.layout, named.semantic_facts.source))
    destination = Agate.Processes.process_id(named) === :invalid_transfer ?
                  :not_a_realized_tracer :
                  only(component_tracers(context.layout, named.semantic_facts.destination))
    parameters = process_parameter_operands(named, context)
    rate = Agate.Compilation.RateOp(
        ExtensionTransferRate(), (input_operand(context.layout, source), parameters.rate)
    )
    return (
        Agate.Compilation.FluxSpec(source, rate, Agate.Compilation.Weight{-1}()),
        Agate.Compilation.FluxSpec(destination, rate, Agate.Compilation.Weight{1}()),
    )
end

@testset "Custom process extension" begin
    components = (source=Agate.Configuration.Pool(:carbon), sink=Agate.Configuration.Pool(:carbon))
    parameters = (transfer_rate=Parameter(0.5),)
    transfer = ExtensionTransfer(:source, :sink; bindings=(rate=:transfer_rate,))

    definition = ModelDefinition(; components, processes=(transfer=transfer,), parameters)
    bgc = Agate.Construction.construct(definition)
    @test bgc(Val(:source), 0, 0, 0, 0, 2.0, 1.0) == -1.0
    @test bgc(Val(:sink), 0, 0, 0, 0, 2.0, 1.0) == 1.0

    invalid_definition = ModelDefinition(;
        components, processes=(invalid_transfer=transfer,), parameters
    )
    message = argument_error_message(() -> Agate.Construction.construct(invalid_definition))
    @test occursin("unrealized targets", message)
    @test occursin(":not_a_realized_tracer", message)
end
