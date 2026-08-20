"""Static setup-time operand that reads a named tracer or auxiliary field."""
struct TracerOp{S} end

"""Static setup-time operand that reads one flattened plankton class."""
struct ClassOp{I} end

"""Static setup-time operand that reads one scalar parameter."""
struct ScalarParamOp{S} end

"""Static setup-time operand that reads one vector parameter element."""
struct VecParamOp{S,I} end

"""Static setup-time operand that reads one matrix parameter element."""
struct MatParamOp{S,I,J} end

"""Static setup-time operand that reads one interaction-matrix element."""
struct InteractionParamOp{S,I,J} end

"""Static operand whose value is one minus another operand."""
struct ComplementOp{O}
    operand::O
end

"""Static operand that evaluates a tuple of child operands."""
struct TupleOp{O}
    operands::O
end

@inline operand_value(::TracerOp{S}, bgc, args) where {S} = getproperty(bgc.tracers, S)(args)
@inline operand_value(::ClassOp{I}, bgc, args) where {I} = bgc.tracers.plankton(args, I)
@inline operand_value(::ScalarParamOp{S}, bgc, args) where {S} = getproperty(bgc.parameters, S)
@inline operand_value(::VecParamOp{S,I}, bgc, args) where {S,I} =
    @inbounds getproperty(bgc.parameters, S)[I]
@inline operand_value(::MatParamOp{S,I,J}, bgc, args) where {S,I,J} =
    @inbounds getproperty(bgc.parameters, S)[I, J]
@inline operand_value(::InteractionParamOp{S,I,J}, bgc, args) where {S,I,J} =
    @inbounds getproperty(bgc.parameters.interactions, S)[I, J]
@inline function operand_value(op::ComplementOp, bgc, args)
    value = operand_value(op.operand, bgc, args)
    return one(value) - value
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

@inline weight_sign(::Weight{Sign}) where {Sign} = Sign

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

"""Setup-time description of one process flux into one concrete tracer."""
struct FluxSpec{R,W}
    process::Symbol
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

function parameter_operand(binding::ParameterBinding, indices::Int...)
    shape = binding.requirement.shape
    path = binding.runtime_path

    if shape === :scalar
        isempty(indices) || throw(ArgumentError("scalar parameter operands do not take indices"))
        length(path) == 1 || throw(
            ArgumentError("scalar runtime parameter paths must have one component; got $path"),
        )
        return ScalarParamOp{only(path)}()
    elseif shape === :vector
        length(indices) == 1 || throw(ArgumentError("vector parameter operands require one index"))
        length(path) == 1 || throw(
            ArgumentError("vector runtime parameter paths must have one component; got $path"),
        )
        return VecParamOp{only(path),only(indices)}()
    elseif shape === :matrix
        length(indices) == 2 || throw(ArgumentError("matrix parameter operands require two indices"))
        i, j = indices
        if length(path) == 1
            return MatParamOp{only(path),i,j}()
        elseif length(path) == 2 && first(path) === :interactions
            return InteractionParamOp{last(path),i,j}()
        end
        throw(ArgumentError("unsupported matrix runtime parameter path $path"))
    end
    throw(ArgumentError("unsupported parameter requirement shape $shape"))
end

function _scalar_component_target(layout::ComponentLayout, component::Symbol)
    hasproperty(layout.component_tracers, component) || throw(
        ArgumentError("unknown scalar component :$component"),
    )
    tracers = getproperty(layout.component_tracers, component)
    length(tracers) == 1 || throw(
        ArgumentError("component :$component must realize to one tracer"),
    )
    return only(tracers)
end

function _realize_population_classes(
    named::NamedProcess,
    populations::Tuple,
    layout::ComponentLayout,
    context::CommunityContext,
)
    tracer_values = ()
    index_values = ()

    for population in populations
        hasproperty(layout.component_tracers, population) || throw(
            ArgumentError("process :$(named.id) references unrealized population :$population"),
        )
        tracers = getproperty(layout.component_tracers, population)
        indices = Tuple(findfirst(==(tracer), context.plankton_symbols) for tracer in tracers)
        any(isnothing, indices) && throw(
            ArgumentError(
                "process :$(named.id) population :$population realizes tracers absent from the current runtime community",
            ),
        )
        tracer_values = (tracer_values..., tracers...)
        index_values = (index_values..., Int.(indices)...)
    end

    return tracer_values, index_values
end

function _target_order(fluxes::Tuple)
    targets = Symbol[]
    for flux in fluxes
        target = flux_target(flux)
        target in targets || push!(targets, target)
    end
    return Tuple(targets)
end

function _validate_target_order(fluxes::Tuple, target_order::Tuple)
    length(unique(target_order)) == length(target_order) || throw(
        ArgumentError("target_order contains duplicate tracer identities"),
    )
    Set(_target_order(fluxes)) == Set(target_order) || throw(
        ArgumentError("target_order must contain exactly the flux targets"),
    )
    return nothing
end

"""Group a flat flux tuple by concrete target tracer."""
function group_fluxes(fluxes::Tuple; target_order=nothing)
    isempty(fluxes) && throw(ArgumentError("cannot group an empty flux tuple"))
    targets = isnothing(target_order) ? _target_order(fluxes) : Tuple(target_order)
    _validate_target_order(fluxes, targets)
    grouped = ntuple(length(targets)) do i
        target = targets[i]
        Tuple(flux for flux in fluxes if flux_target(flux) === target)
    end
    return NamedTuple{targets}(grouped)
end

struct StaticFluxEquation{T}
    terms::T
end

@inline function _sum_fluxes(terms::Tuple{T}, bgc, args) where {T}
    return first(terms)(bgc, args)
end

@inline function _sum_fluxes(terms::Tuple{T,S,Vararg{Any,N}}, bgc, args) where {T,S,N}
    return first(terms)(bgc, args) + _sum_fluxes(Base.tail(terms), bgc, args)
end

@inline function (equation::StaticFluxEquation)(bgc, x, y, z, t, args...)
    return _sum_fluxes(equation.terms, bgc, args)
end

"""Lower one target's setup-time flux tuple to a concrete compiled equation."""
function compile_tendency(fluxes::Tuple)
    isempty(fluxes) && throw(ArgumentError("cannot compile an empty flux tuple"))
    target = flux_target(first(fluxes))
    all(flux -> flux_target(flux) === target, fluxes) || throw(
        ArgumentError("compile_tendency requires fluxes for one target tracer"),
    )
    terms = map(flux -> FluxTerm(flux.rate, flux.weight), fluxes)
    return CompiledEquation(StaticFluxEquation(terms))
end

"""Compile each grouped target flux tuple into a concrete tracer equation."""
compile_tendencies(grouped::NamedTuple) = map(compile_tendency, grouped)

"""Derive all generic process fluxes for a normalized model."""
function model_fluxes(
    definition::NormalizedModelDefinition,
    layout::ComponentLayout,
    context::CommunityContext,
)
    fluxes = ()
    for named in values(definition.processes)
        fluxes = (fluxes..., process_fluxes(named, definition, layout, context)...)
    end
    return fluxes
end

"""Compile a normalized model into one static equation per requested concrete tracer."""
function compile_model_tendencies(
    definition::NormalizedModelDefinition,
    layout::ComponentLayout,
    context::CommunityContext;
    target_order::Tuple,
)
    fluxes = model_fluxes(definition, layout, context)
    grouped = group_fluxes(fluxes; target_order)
    return compile_tendencies(grouped)
end
