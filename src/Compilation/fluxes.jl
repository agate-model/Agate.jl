"""Static setup-time operand that reads a named tracer or auxiliary field."""
struct TracerOp{S} end

"""Static setup-time operand that reads one scalar parameter."""
struct ScalarParamOp{S} end

"""Static setup-time operand that reads one vector parameter element."""
struct VecParamOp{S,I} end

"""Static setup-time operand that reads one matrix parameter element."""
struct MatParamOp{S,I,J} end

"""Static operand whose value is one minus another operand."""
struct ComplementOp{O}
    operand::O
end

"""Static operand whose value is one minus the sum of child operands."""
struct OneMinusSumOp{O}
    operands::O
end

"""Static operand that divides one resolved operand value by another."""
struct RatioOp{N,D}
    numerator::N
    denominator::D
end

"""Static operand that evaluates a tuple of child operands."""
struct TupleOp{O}
    operands::O
end

@inline operand_value(::TracerOp{S}, bgc, args) where {S} = getproperty(bgc.tracers, S)(args)
@inline operand_value(::ScalarParamOp{S}, bgc, args) where {S} = getproperty(bgc.parameters, S)
@inline operand_value(::VecParamOp{S,I}, bgc, args) where {S,I} =
    @inbounds getproperty(bgc.parameters, S)[I]
@inline operand_value(::MatParamOp{S,I,J}, bgc, args) where {S,I,J} =
    @inbounds getproperty(bgc.parameters, S)[I, J]
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
@inline operand_value(op::RatioOp, bgc, args) =
    operand_value(op.numerator, bgc, args) / operand_value(op.denominator, bgc, args)
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
    rank = length(binding.axes)
    parameter = binding.parameter

    if rank == 0
        isempty(indices) || throw(ArgumentError("scalar parameter operands do not take indices"))
        return ScalarParamOp{parameter}()
    elseif rank == 1
        length(indices) == 1 || throw(ArgumentError("vector parameter operands require one index"))
        return VecParamOp{parameter,only(indices)}()
    elseif rank == 2
        length(indices) == 2 || throw(ArgumentError("matrix parameter operands require two indices"))
        i, j = indices
        return MatParamOp{parameter,i,j}()
    end
    throw(ArgumentError("unsupported parameter rank $rank"))
end

@inline _axis_position(local_index::Int, plankton_index::Union{Nothing,Int}=nothing) =
    (; local_index, plankton_index)

function _explicit_storage_index(
    storage_axis::Symbol, position, layout::ModelLayout, parameter::Symbol
)
    plankton_index = position.plankton_index
    isnothing(plankton_index) && throw(
        ArgumentError(
            "parameter :$parameter uses explicit storage axis :$storage_axis for a non-plankton process axis",
        ),
    )
    storage_axis === :plankton && return plankton_index

    storage_indices = axis_indices(layout, storage_axis)
    storage_index = findfirst(==(plankton_index), storage_indices)
    isnothing(storage_index) && throw(
        ArgumentError(
            "parameter :$parameter process class index $plankton_index is not present on storage axis :$storage_axis",
        ),
    )
    return storage_index
end

"""Resolve process-local/global axis positions onto one parameter's runtime storage."""
function parameter_operand(
    binding::ParameterBinding, layout::ModelLayout, axis_positions::NamedTuple
)
    rank = length(binding.axes)
    rank == 0 && return parameter_operand(binding)

    positions = Tuple(
        begin
            hasproperty(axis_positions, axis) || throw(
                ArgumentError(
                    "parameter :$(binding.parameter) axis :$axis has no realized runtime position",
                ),
            )
            getproperty(axis_positions, axis)
        end for axis in binding.axes
    )

    storage_axes = binding.storage_axes
    if storage_axes === nothing
        return parameter_operand(binding, Tuple(position.local_index for position in positions)...)
    elseif rank == 1
        storage_axes isa Symbol || throw(
            ArgumentError(
                "vector parameter :$(binding.parameter) must have one Symbol storage axis",
            ),
        )
        index = _explicit_storage_index(storage_axes, only(positions), layout, binding.parameter)
        return parameter_operand(binding, index)
    elseif rank == 2
        storage_axes isa Tuple && length(storage_axes) == 2 || throw(
            ArgumentError(
                "matrix parameter :$(binding.parameter) must have two storage axes",
            ),
        )
        indices = ntuple(2) do i
            _explicit_storage_index(storage_axes[i], positions[i], layout, binding.parameter)
        end
        return parameter_operand(binding, indices...)
    end

    throw(ArgumentError("unsupported parameter rank $rank"))
end

function _scalar_component_target(layout::ModelLayout, component::Symbol)
    hasproperty(layout.component_tracers, component) || throw(
        ArgumentError("unknown scalar component :$component"),
    )
    tracers = getproperty(layout.component_tracers, component)
    length(tracers) == 1 || throw(
        ArgumentError("component :$component must realize to one tracer"),
    )
    return only(tracers)
end

"""Resolve one explicit population state onto static tracer operands for its ecological classes."""
function _realize_population_state(
    named::NamedProcess,
    reference::PopulationStateRef,
    layout::ModelLayout,
)
    population = reference.population
    hasproperty(layout.component_classes, population) || throw(
        ArgumentError("process :$(named.id) references unrealized population :$population"),
    )
    classes = component_classes(layout, population)
    tracers = Tuple(state_tracer(layout, reference, i) for i in eachindex(classes))
    return tracers, component_class_indices(layout, population)
end

"""Resolve one population state read at a global ecological class index to a static tracer operand."""
function state_operand(
    layout::ModelLayout,
    reference::PopulationStateRef,
    global_class_index::Integer,
)
    1 <= global_class_index <= length(layout.class_symbols) || throw(
        ArgumentError("global ecological class index $global_class_index is out of bounds"),
    )
    index = Int(global_class_index)
    layout.class_populations[index] === reference.population || throw(
        ArgumentError(
            "ecological class :$(layout.class_symbols[index]) does not belong to population :$(reference.population)",
        ),
    )
    tracer = state_tracer(layout, reference, layout.component_local_indices[index])
    return TracerOp{tracer}()
end

function _realize_population_classes(
    named::NamedProcess,
    populations::Tuple,
    layout::ModelLayout,
)
    tracer_values = Symbol[]
    index_values = Int[]

    for population in populations
        state_mapping = component_state_tracers(layout, population)
        state_mapping isa NamedTuple || throw(
            ArgumentError("process :$(named.id) component :$population is not a Population"),
        )
        length(state_mapping) == 1 || throw(
            ArgumentError(
                "process :$(named.id) requires explicit state selection for multi-state population :$population",
            ),
        )
        reference = PopulationStateRef(population, only(keys(state_mapping)))
        tracers, indices = _realize_population_state(named, reference, layout)
        append!(tracer_values, tracers)
        append!(index_values, indices)
    end

    return Tuple(tracer_values), Tuple(index_values)
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
    issubset(Set(_target_order(fluxes)), Set(target_order)) || throw(
        ArgumentError("target_order must contain every flux target"),
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

"""Lower one target's setup-time flux tuple to a concrete compiled equation."""
function compile_tendency(fluxes::Tuple)
    isempty(fluxes) && return CompiledEquation(StaticZeroEquation())
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
    layout::ModelLayout,
)
    fluxes = Any[]
    for named in values(definition.processes)
        append!(fluxes, process_fluxes(named, definition, layout))
    end
    return Tuple(fluxes)
end

"""Compile a normalized model into one static equation per requested concrete tracer."""
function compile_model_tendencies(
    definition::NormalizedModelDefinition,
    layout::ModelLayout;
    target_order::Tuple,
)
    fluxes = model_fluxes(definition, layout)
    grouped = group_fluxes(fluxes; target_order)
    return compile_tendencies(grouped)
end
