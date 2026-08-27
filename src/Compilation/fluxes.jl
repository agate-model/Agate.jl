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

@inline _axis_position(local_index::Int, class::Symbol) = (; local_index, class)

"""Resolve one tracer or auxiliary identity to its final positional runtime input."""
function input_operand(layout::ModelLayout, identity::Symbol)
    hasproperty(layout.input_indices, identity) || throw(
        ArgumentError("unknown realized tracer or auxiliary input :$identity"),
    )
    return InputOp{getproperty(layout.input_indices, identity)}()
end

"""Resolve one canonical parameter slot directly from realized ecological class identity."""
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
        parameter_storage_index(plan, binding.parameter, dimension, position.class)
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
    refs isa NamedTuple || throw(ArgumentError(
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

function _realize_population_states(references::Tuple, layout::ModelLayout)
    tracers = Symbol[]
    for reference in references
        classes = component_classes(layout, reference.population)
        append!(tracers, (state_tracer(layout, reference, i) for i in eachindex(classes)))
    end
    return Tuple(tracers)
end

function _realize_population_classes(references::Tuple, layout::ModelLayout)
    classes = Symbol[]
    for reference in references
        append!(classes, component_classes(layout, reference.population))
    end
    return Tuple(classes)
end

function _realize_component_classes(components::Tuple, layout::ModelLayout)
    classes = Symbol[]
    for component in components
        append!(classes, component_classes(layout, component))
    end
    return Tuple(classes)
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
