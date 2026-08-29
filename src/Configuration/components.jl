"""Population state described by intrinsic component properties.

Ecological function is supplied by process participation. `states` names every prognostic
state carried by the population, while `reference_state` identifies the state used as its
biological reference basis. Elemental bookkeeping is inferred centrally by [`state_element`](@ref);
non-elemental states such as `:chlorophyll` return `nothing`.

A one-state population may be authored with its reference state as the positional argument.
"""
struct Population{ST<:Tuple,S}
    states::ST
    reference_state::Symbol
    size_structure::S

    function Population(states::ST, reference_state::Symbol, size_structure::S) where {ST<:Tuple,S}
        isempty(states) && throw(ArgumentError("Population must define at least one state."))
        all(state -> state isa Symbol, states) || throw(
            ArgumentError("Population `states` must contain only Symbols."),
        )
        length(unique(states)) == length(states) || throw(
            ArgumentError("Population state identities must be unique."),
        )
        reference_state in states || throw(
            ArgumentError("Population `reference_state` must be one of its declared states."),
        )
        return new{ST,S}(states, reference_state, size_structure)
    end
end

function _canonical_states(states)
    states isa Symbol && return (states,)
    states isa Tuple || throw(
        ArgumentError("Population `states` must be a Symbol or Tuple of Symbols."),
    )
    return states
end

function Population(reference_state; size_structure=nothing)
    reference_state isa Symbol || throw(
        ArgumentError("Population reference state must be a Symbol."),
    )
    return Population((reference_state,), reference_state, size_structure)
end

function Population(; states, reference_state, size_structure=nothing)
    reference_state isa Symbol || throw(
        ArgumentError("Population `reference_state` must be a Symbol."),
    )
    return Population(_canonical_states(states), reference_state, size_structure)
end

"""Material-pool state described by intrinsic component properties."""
struct Pool{C,S}
    currency::C
    size_structure::S
end

function Pool(currency; size_structure=nothing)
    isnothing(currency) && throw(ArgumentError("Pool currency must be specified."))
    return Pool(currency, size_structure)
end

"""Reference one named prognostic state carried by a logical population."""
struct PopulationStateRef
    population::Symbol
    state::Symbol
end

"""Return all prognostic state identities for a population."""
@inline states(component::Population) = component.states

"""Return the population state used as its biological reference basis."""
@inline reference_state(component::Population) = component.reference_state

"""Return all prognostic population states other than the reference state."""
@inline variable_states(component::Population) = Tuple(
    state for state in component.states if state !== component.reference_state
)

const _ELEMENTAL_STATES = (:carbon, :nitrogen, :phosphorus, :iron, :silicon)

"""Return the conserved Element represented by `state`, or `nothing` for a non-elemental state.

Element identity is the additive elemental inventory represented by the state's numerical value,
not the chemical composition of the represented material. Canonical elemental state names map to
themselves; other states remain non-elemental until an explicit mapping API is introduced.
"""
function state_element(component::Population, state::Symbol)
    state in component.states || throw(ArgumentError("Population has no state :$state."))
    return state in _ELEMENTAL_STATES ? state : nothing
end

"""Return the conserved-material currency represented by a material Pool."""
@inline currency(component::Pool) = component.currency

"""Return the intrinsic size-structure specification for `component`."""
@inline size_structure(component::Union{Population,Pool}) = component.size_structure
