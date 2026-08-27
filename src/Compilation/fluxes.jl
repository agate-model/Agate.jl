"""Setup-time context shared by process lowering and static tendency compilation."""
struct CompileContext{D<:NormalizedModelDefinition,L<:ModelLayout,P<:ParameterPlan}
    definition::D
    layout::L
    plan::P
end

"""Lower one normalized process to generic fluxes.

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

@inline _axis_position(local_index::Int) = (; local_index)

"""Resolve one tracer or auxiliary identity to its final positional runtime input."""
function input_operand(layout::ModelLayout, identity::Symbol)
    hasproperty(layout.input_indices, identity) || throw(
        ArgumentError("unknown realized tracer or auxiliary input :$identity"),
    )
    return InputOp{getproperty(layout.input_indices, identity)}()
end

"""Resolve one normalized parameter slot through its precomputed `ParameterPlan` mapping."""
function parameter_operand(
    binding::ParameterBinding,
    plan::ParameterPlan,
    axis_positions::NamedTuple=NamedTuple(),
)
    rank = length(binding.axes)
    indices = if rank == 0
        ()
    else
        slot = planned_parameter_slot(plan, binding)
        ntuple(rank) do dimension
            axis = binding.axes[dimension]
            hasproperty(axis_positions, axis) || throw(ArgumentError(
                "parameter :$(binding.parameter) axis :$axis has no realized runtime position",
            ))
            local_index = getproperty(axis_positions, axis).local_index
            mapping = slot[dimension]
            1 <= local_index <= length(mapping) || throw(ArgumentError(
                "parameter :$(binding.parameter) axis :$axis local index $local_index is out of bounds",
            ))
            mapping[local_index]
        end
    end
    return ParameterOp{binding.parameter,indices}()
end

function _scalar_component_target(layout::ModelLayout, component::Symbol)
    tracers = getproperty(layout.component_tracers, component)
    return only(tracers)
end

function _realize_normalized_population_states(references::Tuple, layout::ModelLayout)
    tracers = Symbol[]
    for reference in references
        classes = component_classes(layout, reference.population)
        append!(tracers, (state_tracer(layout, reference, i) for i in eachindex(classes)))
    end
    return Tuple(tracers)
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

"""Derive all generic process fluxes for a normalized model."""
function model_fluxes(context::CompileContext)
    fluxes = Any[]
    for named in values(context.definition.processes)
        append!(fluxes, process_fluxes(named, context))
    end
    return Tuple(fluxes)
end

"""Compile a normalized model into one static equation per requested concrete tracer."""
function compile_model_tendencies(context::CompileContext; target_order::Tuple)
    return compile_tendencies(group_fluxes(model_fluxes(context); target_order))
end
