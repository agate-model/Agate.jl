using ForwardDiff
using Test

using Agate.Configuration:
    Population, Pool, population_state, state_tracer
using Agate.Construction: construct
using Agate.Parameters:
    ParameterDefinition, ParameterProvision, FillDefault
using Agate.Processes:
    AbstractProcess, AbstractFormulation, ModelDefinition, ParameterSlot,
    formulation, parameter_slot_bindings, process_id
using Agate.Compilation:
    TracerOp, RatioOp, InteractionParamOp, RateElement, Weight, FluxSpec,
    parameter_operand, state_operand, _axis_position, _realize_population_state

import Agate.Processes:
    formulation_tag, parameter_slots, process_kind, participants, rate_axes,
    process_rate, uses_living_interactions
import Agate.Compilation: process_fluxes

struct LinearStateTurnover <: AbstractFormulation end

struct StateTurnover{D<:NamedTuple} <: AbstractProcess
    formulation::LinearStateTurnover
    population::Symbol
    reference_state::Symbol
    destinations::D
end

StateTurnover(population::Symbol, reference_state::Symbol, destinations::NamedTuple) =
    StateTurnover(LinearStateTurnover(), population, reference_state, destinations)

formulation_tag(::LinearStateTurnover) = :linear_state_turnover
parameter_slots(::LinearStateTurnover) = (ParameterSlot(:rate, (:population,)),)
process_kind(::StateTurnover) = :state_turnover
participants(process::StateTurnover) = (
    population=(process.population,), destination=Tuple(values(process.destinations))
)
rate_axes(::StateTurnover) = (:population,)
process_rate(::LinearStateTurnover, reference_state, rate) = reference_state * rate

function process_fluxes(
    named::Agate.Processes.NamedProcess{P},
    definition::Agate.Processes.NormalizedModelDefinition,
    layout::Agate.Configuration.ComponentLayout,
    context::Agate.Configuration.CommunityContext,
) where {P<:StateTurnover}
    process = named.process
    reference = population_state(process.population, process.reference_state)
    reference_tracers, population_indices = _realize_population_state(
        named, reference, layout, context
    )
    slots = parameter_slot_bindings(definition, named, (), formulation(process))
    fluxes = ()

    for population_axis in eachindex(reference_tracers)
        population_index = population_indices[population_axis]
        axis_positions = (
            population=_axis_position(population_axis, population_index),
        )
        reference_operand = TracerOp{reference_tracers[population_axis]}()
        rate = RateElement(
            formulation(process),
            (
                reference_operand,
                parameter_operand(slots.rate, context, axis_positions),
            ),
        )

        for (state, destination) in pairs(process.destinations)
            state_reference = population_state(process.population, state)
            source = state_tracer(layout, state_reference, population_axis)
            source_operand = state_operand(
                layout, context, state_reference, population_index
            )
            ratio = state === process.reference_state ? () :
                    (RatioOp(source_operand, reference_operand),)
            fluxes = (
                fluxes...,
                FluxSpec(process_id(named), source, rate, Weight{-1}(ratio)),
                FluxSpec(process_id(named), destination, rate, Weight{1}(ratio)),
            )
        end
    end
    return fluxes
end

struct StateInteractionRate <: AbstractFormulation end

struct StateInteraction <: AbstractProcess
    formulation::StateInteractionRate
    consumer::Symbol
    resource::Symbol
    carbon_state::Symbol
    nitrogen_state::Symbol
end

StateInteraction(consumer, resource; carbon_state=:carbon, nitrogen_state=:nitrogen) =
    StateInteraction(StateInteractionRate(), consumer, resource, carbon_state, nitrogen_state)

formulation_tag(::StateInteractionRate) = :state_interaction
parameter_slots(::StateInteractionRate) = (
    ParameterSlot(:palatability, (:consumer, :resource)),
)
uses_living_interactions(::StateInteractionRate) = true
process_kind(::StateInteraction) = :state_interaction
participants(process::StateInteraction) = (
    consumer=(process.consumer,), resource=(process.resource,)
)
rate_axes(::StateInteraction) = (:consumer, :resource)
process_rate(::StateInteractionRate, resource, consumer, palatability) =
    palatability * resource * consumer

function process_fluxes(
    named::Agate.Processes.NamedProcess{P},
    definition::Agate.Processes.NormalizedModelDefinition,
    layout::Agate.Configuration.ComponentLayout,
    context::Agate.Configuration.CommunityContext,
) where {P<:StateInteraction}
    process = named.process
    consumer_carbon = population_state(process.consumer, process.carbon_state)
    resource_carbon = population_state(process.resource, process.carbon_state)
    consumer_tracers, consumer_indices = _realize_population_state(
        named, consumer_carbon, layout, context
    )
    resource_tracers, resource_indices = _realize_population_state(
        named, resource_carbon, layout, context
    )
    slots = parameter_slot_bindings(definition, named, (), formulation(process))
    fluxes = ()

    for consumer_axis in eachindex(consumer_tracers)
        consumer_index = consumer_indices[consumer_axis]
        consumer_carbon_operand = TracerOp{consumer_tracers[consumer_axis]}()
        for resource_axis in eachindex(resource_tracers)
            resource_index = resource_indices[resource_axis]
            resource_carbon_operand = TracerOp{resource_tracers[resource_axis]}()
            axis_positions = (
                consumer=_axis_position(consumer_axis, consumer_index),
                resource=_axis_position(resource_axis, resource_index),
            )
            palatability = parameter_operand(slots.palatability, context, axis_positions)
            rate = RateElement(
                formulation(process),
                (resource_carbon_operand, consumer_carbon_operand, palatability),
            )

            resource_nitrogen = population_state(process.resource, process.nitrogen_state)
            consumer_nitrogen = population_state(process.consumer, process.nitrogen_state)
            resource_nitrogen_operand = state_operand(
                layout, context, resource_nitrogen, resource_index
            )
            nitrogen_ratio = RatioOp(
                resource_nitrogen_operand, resource_carbon_operand
            )
            resource_nitrogen_tracer = state_tracer(
                layout, resource_nitrogen, resource_axis
            )
            consumer_nitrogen_tracer = state_tracer(
                layout, consumer_nitrogen, consumer_axis
            )

            fluxes = (
                fluxes...,
                FluxSpec(
                    process_id(named), resource_tracers[resource_axis], rate, Weight{-1}()
                ),
                FluxSpec(
                    process_id(named), consumer_tracers[consumer_axis], rate, Weight{1}()
                ),
                FluxSpec(
                    process_id(named), resource_nitrogen_tracer, rate,
                    Weight{-1}((nitrogen_ratio,)),
                ),
                FluxSpec(
                    process_id(named), consumer_nitrogen_tracer, rate,
                    Weight{1}((nitrogen_ratio,)),
                ),
            )
        end
    end
    return fluxes
end

function multistate_args(bgc, values::NamedTuple)
    names = Agate.Introspection.tracer_names(bgc)
    return (0.0, 0.0, 0.0, 0.0, Tuple(getproperty(values, name) for name in names)...)
end

@testset "Multi-state state-aware compilation and conservation" begin
    components = (
        DOC=Pool(:carbon),
        DON=Pool(:nitrogen),
        DOP=Pool(:phosphorus),
        P=Population(;
            states=(carbon=:carbon, nitrogen=:nitrogen, phosphorus=:phosphorus),
            size_structure=[1.0, 2.0],
        ),
    )
    process = StateTurnover(
        :P, :carbon, (carbon=:DOC, nitrogen=:DON, phosphorus=:DOP)
    )
    parameters = (
        ParameterDefinition(
            :turnover_rate,
            FillDefault(0.1);
            shape=:vector,
            axes=:plankton,
            provides=ParameterProvision(:turnover, :rate),
        ),
    )
    bgc = construct(ModelDefinition(;
        components, processes=(turnover=process,), parameters
    ))

    @test Agate.Introspection.tracer_names(bgc) == [
        :DOC, :DON, :DOP,
        :P_1_carbon, :P_1_nitrogen, :P_1_phosphorus,
        :P_2_carbon, :P_2_nitrogen, :P_2_phosphorus,
    ]
    @test bgc.plankton_diameters == (1.0, 2.0)
    @test bgc.tracers.idx.plankton_base == 0

    values = (
        DOC=0.0, DON=0.0, DOP=0.0,
        P_1_carbon=10.0, P_1_nitrogen=2.0, P_1_phosphorus=1.0,
        P_2_carbon=20.0, P_2_nitrogen=6.0, P_2_phosphorus=4.0,
    )
    args = multistate_args(bgc, values)
    tendencies = NamedTuple{Tuple(Agate.Introspection.tracer_names(bgc))}(
        Tuple(bgc(Val(name), args...) for name in Agate.Introspection.tracer_names(bgc))
    )

    @test tendencies.P_1_carbon ≈ -1.0
    @test tendencies.P_1_nitrogen ≈ -0.2
    @test tendencies.P_1_phosphorus ≈ -0.1
    @test tendencies.P_2_carbon ≈ -2.0
    @test tendencies.P_2_nitrogen ≈ -0.6
    @test tendencies.P_2_phosphorus ≈ -0.4
    @test tendencies.DOC ≈ 3.0
    @test tendencies.DON ≈ 0.8
    @test tendencies.DOP ≈ 0.5
    @test tendencies.DOC + tendencies.P_1_carbon + tendencies.P_2_carbon ≈ 0.0
    @test tendencies.DON + tendencies.P_1_nitrogen + tendencies.P_2_nitrogen ≈ 0.0
    @test tendencies.DOP + tendencies.P_1_phosphorus + tendencies.P_2_phosphorus ≈ 0.0
    @test all(equation -> isbitstype(typeof(equation)), values(bgc.tracer_functions))

    derivative = ForwardDiff.derivative(2.0) do nitrogen
        dynamic = merge(values, (P_1_nitrogen=nitrogen,))
        bgc(Val(:P_1_nitrogen), multistate_args(bgc, dynamic)...)
    end
    @test derivative ≈ -0.1
end

@testset "One ecological interaction spans multiple population states" begin
    components = (
        Z=Population(; states=(carbon=:carbon, nitrogen=:nitrogen)),
        P=Population(; states=(carbon=:carbon, nitrogen=:nitrogen)),
    )
    process = StateInteraction(:Z, :P)
    parameters = (
        ParameterDefinition(
            :state_palatability,
            FillDefault(0.5);
            shape=:matrix,
            axes=(:consumer, :prey),
            runtime_path=(:interactions, :state_palatability),
            provides=ParameterProvision(:consume, :palatability),
        ),
    )
    bgc = construct(ModelDefinition(;
        components, processes=(consume=process,), parameters
    ))

    @test size(bgc.parameters.interactions.state_palatability) == (1, 1)
    @test bgc.parameters.interactions.state_palatability[1, 1] == 0.5

    values = (
        Z_carbon=1.0, Z_nitrogen=0.2,
        P_carbon=2.0, P_nitrogen=1.0,
    )
    args = multistate_args(bgc, values)
    @test bgc(Val(:P_carbon), args...) ≈ -1.0
    @test bgc(Val(:Z_carbon), args...) ≈ 1.0
    @test bgc(Val(:P_nitrogen), args...) ≈ -0.5
    @test bgc(Val(:Z_nitrogen), args...) ≈ 0.5
    @test bgc(Val(:P_carbon), args...) + bgc(Val(:Z_carbon), args...) ≈ 0.0
    @test bgc(Val(:P_nitrogen), args...) + bgc(Val(:Z_nitrogen), args...) ≈ 0.0

    normalized = Agate.Processes.normalize_model(ModelDefinition(;
        components, processes=(consume=process,), parameters
    ))
    realization = Agate.Construction._realize_process_definition(normalized, Float64)
    fluxes = process_fluxes(
        normalized.processes.consume, normalized, realization.layout, realization.context
    )
    interaction_operand = fluxes[1].rate.operands[3]
    @test interaction_operand == InteractionParamOp{:state_palatability,1,1}()
    @test fluxes[3].weight.operands[1] isa RatioOp
end
