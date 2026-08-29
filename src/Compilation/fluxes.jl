"""Setup-time context shared by process lowering and static tendency compilation."""
struct CompileContext{D<:CanonicalModelDefinition,L<:ModelLayout,P<:ParameterPlan}
    definition::D
    layout::L
    plan::P
end

"""Lower one canonical process to generic fluxes.

Custom process implementations may extend this hook for a concrete `NamedProcess` subtype.
"""
function process_fluxes(named::NamedProcess, context::CompileContext)
    throw(ArgumentError(
        "no compiler lowering is defined for process type $(typeof(named.process))",
    ))
end

"""Setup-time description of one flux into one concrete tracer."""
struct FluxSpec{R,W}
    target::Symbol
    rate::R
    weight::W
end

@inline flux_target(flux::FluxSpec) = flux.target

@inline _axis_position(axis_index::Int, entity::Symbol, component_index::Int) =
    (; axis_index, entity, component_index)

"""Realize participant states or components to concrete tracers and ecological positions."""
function _realize_participants(items::Tuple, layout::ModelLayout)
    participants = Any[]
    axis_index = 0
    for item in items
        is_state_reference = item isa PlanktonStateRef
        component = is_state_reference ? item.plankton : item
        tracers = is_state_reference ? state_tracers(layout, item) : component_tracers(layout, item)
        entities = component_entities(layout, component)
        length(tracers) == length(entities) || throw(ArgumentError(
            "participant $item must realize exactly one tracer per realized entity",
        ))
        for component_index in eachindex(entities)
            axis_index += 1
            position = _axis_position(axis_index, entities[component_index], component_index)
            push!(participants, (; tracer=tracers[component_index], component, component_index, position))
        end
    end
    return Tuple(participants)
end

"""Resolve one tracer or auxiliary identity to its final positional runtime input."""
function input_operand(layout::ModelLayout, identity::Symbol)
    hasproperty(layout.input_indices, identity) || throw(
        ArgumentError("unknown realized tracer or auxiliary input :$identity"),
    )
    return InputOp{getproperty(layout.input_indices, identity)}()
end

"""Resolve one canonical parameter slot directly from realized entity identity."""
function parameter_operand(
    binding::ParameterBinding,
    plan::ParameterPlan,
    axis_positions::NamedTuple=NamedTuple(),
)
    rank = length(binding.axes)
    indices = ntuple(rank) do dimension
        axis = binding.axes[dimension]
        hasproperty(axis_positions, axis) || throw(ArgumentError(
            "parameter :$(binding.parameter) axis :$axis has no realized runtime position",
        ))
        position = getproperty(axis_positions, axis)
        parameter_storage_index(plan, binding.parameter, dimension, position.entity)
    end
    return ParameterOp{binding.parameter,indices}()
end


"""Resolve one dense canonical binding reference to its static runtime parameter operand."""
function parameter_operand(
    ref::Int,
    context::CompileContext,
    axis_positions::NamedTuple=NamedTuple(),
)
    binding = context.definition.parameter_bindings[ref]
    return parameter_operand(binding, context.plan, axis_positions)
end

"""Return compiled operands for the unqualified process-owned parameter slots of a custom process.

This is the narrow parameterized-process extension seam: custom lowering receives semantic slot
names mapped directly to static operands without depending on canonical binding references.
"""
function process_parameter_operands(
    named::NamedProcess,
    context::CompileContext,
    axis_positions::NamedTuple=NamedTuple(),
)
    refs = named.binding_refs.process
    all(ref -> ref isa Int, values(refs)) || throw(ArgumentError(
        "process :$(process_id(named)) has qualifier-specific parameter slots; custom lowering must resolve its setup facts explicitly",
    ))
    names = keys(refs)
    return NamedTuple{names}(Tuple(
        parameter_operand(ref, context, axis_positions) for ref in values(refs)
    ))
end

function _scalar_component_target(layout::ModelLayout, component::Symbol)
    tracers = getproperty(layout.component_tracers, component)
    return only(tracers)
end

function _target_order(fluxes::Tuple)
    targets = Symbol[]
    for flux in fluxes
        target = flux_target(flux)
        target in targets || push!(targets, target)
    end
    return Tuple(targets)
end

"""Group a flat flux tuple by concrete target tracer."""
function group_fluxes(fluxes::Tuple; target_order=nothing)
    targets = isnothing(target_order) ? _target_order(fluxes) : Tuple(target_order)
    if !isnothing(target_order)
        unknown = Tuple(target for target in _target_order(fluxes) if !(target in targets))
        isempty(unknown) || throw(ArgumentError(
            "process lowering produced fluxes for unrealized targets $unknown; " *
            "targets must be concrete tracers in the realized layout",
        ))
    end
    grouped = ntuple(length(targets)) do i
        target = targets[i]
        Tuple(flux for flux in fluxes if flux_target(flux) === target)
    end
    return NamedTuple{targets}(grouped)
end

"""Lower one target's grouped setup-time flux tuple to a concrete runtime equation."""
function compile_tendency(fluxes::Tuple)
    isempty(fluxes) && return StaticZeroEquation()
    terms = map(flux -> FluxTerm(flux.rate, flux.weight), fluxes)
    return StaticFluxEquation(terms)
end

"""Compile each grouped target flux tuple into a concrete tracer equation."""
compile_tendencies(grouped::NamedTuple) = map(compile_tendency, grouped)

"""Derive all generic process fluxes for a canonical model."""
function model_fluxes(context::CompileContext)
    fluxes = Any[]
    for named in values(context.definition.processes)
        append!(fluxes, process_fluxes(named, context))
    end
    return Tuple(fluxes)
end

"""Compile a canonical model into one static equation per requested concrete tracer."""
function compile_model_tendencies(context::CompileContext; target_order::Tuple)
    return compile_tendencies(group_fluxes(model_fluxes(context); target_order))
end
