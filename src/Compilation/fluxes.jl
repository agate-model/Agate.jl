"""Static operand that reads one pre-indexed tracer or auxiliary input."""
struct InputOp{S,I} end

"""Static operand that reads one pre-indexed runtime parameter value."""
struct ParameterOp{S,I} end

"""Static operand whose value is one minus another operand."""
struct ComplementOp{O}
    operand::O
end

"""Static operand whose value is one minus the sum of child operands."""
struct OneMinusSumOp{O}
    operands::O
end

"""Static operand that evaluates a tuple of child operands."""
struct TupleOp{O}
    operands::O
end

@inline operand_value(::InputOp{S,I}, bgc, args) where {S,I} = @inbounds args[I]
@inline operand_value(::ParameterOp{S,()}, bgc, args) where {S} = getproperty(bgc.parameters, S)
@inline operand_value(::ParameterOp{S,I}, bgc, args) where {S,I} =
    @inbounds getproperty(bgc.parameters, S)[I...]

@inline function operand_value(op::ComplementOp, bgc, args)
    value = operand_value(op.operand, bgc, args)
    return one(value) - value
end

@inline function _operand_sum(operands::Tuple{T}, bgc, args) where {T}
    return operand_value(first(operands), bgc, args)
end

@inline function _operand_sum(operands::Tuple{T,S,Vararg{Any,N}}, bgc, args) where {T,S,N}
    return operand_value(first(operands), bgc, args) + _operand_sum(Base.tail(operands), bgc, args)
end

@inline function operand_value(op::OneMinusSumOp, bgc, args)
    total = _operand_sum(op.operands, bgc, args)
    return one(total) - total
end
@inline operand_value(op::TupleOp, bgc, args) = operand_values(op.operands, bgc, args)

@inline operand_values(::Tuple{}, bgc, args) = ()
@inline function operand_values(operands::Tuple, bgc, args)
    return (operand_value(first(operands), bgc, args), operand_values(Base.tail(operands), bgc, args)...)
end

"""One multiplicative factor attached to a canonical process rate."""
struct FactorElement{F,O}
    formulation::F
    operands::O
end

@inline function factor_element_value(factor::FactorElement, bgc, args)
    values = operand_values(factor.operands, bgc, args)
    return factor_value(factor.formulation, values...)
end

@inline operand_value(factor::FactorElement, bgc, args) =
    factor_element_value(factor, bgc, args)

@inline apply_rate_factors(::Tuple{}, bgc, args, rate) = rate
@inline function apply_rate_factors(factors::Tuple, bgc, args, rate)
    value = factor_element_value(first(factors), bgc, args)
    return apply_rate_factors(Base.tail(factors), bgc, args, rate * value)
end

"""Canonical runtime process rate: formulation plus static operand tuple."""
struct RateElement{F,O,X}
    formulation::F
    operands::O
    factors::X
end

RateElement(formulation, operands::Tuple; factors=()) =
    RateElement{typeof(formulation),typeof(operands),typeof(factors)}(formulation, operands, factors)

@inline function (rate::RateElement)(bgc, args)
    values = operand_values(rate.operands, bgc, args)
    value = process_rate(rate.formulation, values...)
    return apply_rate_factors(rate.factors, bgc, args, value)
end

"""Static multiplicative flux weight with sign encoded in the type."""
struct Weight{Sign,O}
    operands::O
end

Weight{Sign}(operands::Tuple=()) where {Sign} = Weight{Sign,typeof(operands)}(operands)

@inline function weight_product(operands::Tuple{T}, bgc, args) where {T}
    return operand_value(first(operands), bgc, args)
end

@inline function weight_product(operands::Tuple{T,S,Vararg{Any,N}}, bgc, args) where {T,S,N}
    return operand_value(first(operands), bgc, args) * weight_product(Base.tail(operands), bgc, args)
end

@inline apply_weight(::Weight{1,Tuple{}}, bgc, args, rate) = rate
@inline apply_weight(::Weight{-1,Tuple{}}, bgc, args, rate) = -rate
@inline apply_weight(weight::Weight{1}, bgc, args, rate) = weight_product(weight.operands, bgc, args) * rate
@inline apply_weight(weight::Weight{-1}, bgc, args, rate) = -weight_product(weight.operands, bgc, args) * rate

"""Setup-time description of one flux into one concrete tracer."""
struct FluxSpec{R,W}
    target::Symbol
    rate::R
    weight::W
end

@inline flux_target(flux::FluxSpec) = flux.target

"""Lean runtime flux term containing only its scientific rate and weight."""
struct FluxTerm{R,W}
    rate::R
    weight::W
end

@inline function (term::FluxTerm)(bgc, args)
    rate = term.rate(bgc, args)
    return apply_weight(term.weight, bgc, args, rate)
end

@inline _axis_position(local_index::Int) = (; local_index)

"""Resolve one tracer or auxiliary identity to its final positional runtime input."""
function input_operand(layout::ModelLayout, identity::Symbol)
    if hasproperty(layout.tracer_indices, identity)
        return InputOp{identity,getproperty(layout.tracer_indices, identity)}()
    elseif hasproperty(layout.auxiliary_indices, identity)
        return InputOp{identity,getproperty(layout.auxiliary_indices, identity)}()
    end
    throw(ArgumentError("unknown realized tracer or auxiliary input :$identity"))
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

struct StaticFluxEquation{T}
    terms::T
end

struct StaticZeroEquation end

@inline (::StaticZeroEquation)(bgc, x, y, z, t, args...) = zero(first(args))

@inline function _sum_fluxes(terms::Tuple{T}, bgc, args) where {T}
    return first(terms)(bgc, args)
end

@inline function _sum_fluxes(terms::Tuple{T,S,Vararg{Any,N}}, bgc, args) where {T,S,N}
    return first(terms)(bgc, args) + _sum_fluxes(Base.tail(terms), bgc, args)
end

@inline function (equation::StaticFluxEquation)(bgc, x, y, z, t, args...)
    return _sum_fluxes(equation.terms, bgc, args)
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
function model_fluxes(
    definition::NormalizedModelDefinition,
    layout::ModelLayout,
    plan::ParameterPlan,
)
    fluxes = Any[]
    for named in values(definition.processes)
        append!(fluxes, process_fluxes(named, definition, layout, plan))
    end
    return Tuple(fluxes)
end

"""Compile a normalized model into one static equation per requested concrete tracer."""
function compile_model_tendencies(
    definition::NormalizedModelDefinition,
    layout::ModelLayout,
    plan::ParameterPlan;
    target_order::Tuple,
)
    fluxes = model_fluxes(definition, layout, plan)
    return compile_tendencies(group_fluxes(fluxes; target_order))
end
