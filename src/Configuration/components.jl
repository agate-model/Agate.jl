"""Plankton state described by intrinsic component properties.

Ecological function is supplied by process participation. `states` names every prognostic
state carried by the plankton, while `reference_state` identifies the state used as its
biological reference basis. Elemental bookkeeping is inferred centrally by [`state_element`](@ref);
non-elemental states such as `:chlorophyll` return `nothing`.

"""
struct Plankton{ST<:Tuple,S}
    states::ST
    reference_state::Symbol
    size_structure::S

    function Plankton(states::ST, reference_state::Symbol, size_structure::S) where {ST<:Tuple,S}
        isempty(states) && throw(ArgumentError("Plankton must define at least one state."))
        all(state -> state isa Symbol, states) || throw(
            ArgumentError("Plankton `states` must contain only Symbols."),
        )
        length(unique(states)) == length(states) || throw(
            ArgumentError("Plankton state identities must be unique."),
        )
        reference_state in states || throw(
            ArgumentError("Plankton `reference_state` must be one of its declared states."),
        )
        return new{ST,S}(states, reference_state, size_structure)
    end
end

function _canonical_states(states)
    states isa Symbol && return (states,)
    states isa Tuple || throw(
        ArgumentError("Plankton `states` must be a Symbol or Tuple of Symbols."),
    )
    return states
end

function Plankton(; states, reference_state, size_structure=nothing)
    reference_state isa Symbol || throw(
        ArgumentError("Plankton `reference_state` must be a Symbol."),
    )
    return Plankton(_canonical_states(states), reference_state, size_structure)
end

"""Scalar material-pool state described by its conserved Element."""
struct Pool
    element::Symbol
end

"""Reference one named prognostic state carried by a logical plankton."""
struct PlanktonStateRef
    plankton::Symbol
    state::Symbol
end

"""Return all prognostic state identities for a plankton."""
@inline states(component::Plankton) = component.states

"""Return the plankton state used as its biological reference basis."""
@inline reference_state(component::Plankton) = component.reference_state

"""Return all prognostic plankton states other than the reference state."""
@inline variable_states(component::Plankton) = Tuple(
    state for state in component.states if state !== component.reference_state
)

const _ELEMENTAL_STATES = (:carbon, :nitrogen, :phosphorus, :iron, :silicon)

"""Return the conserved Element represented by `state`, or `nothing` for a non-elemental state.

Element identity is the additive elemental inventory represented by the state's numerical value,
not the chemical composition of the represented material. Canonical elemental state names map to
themselves; other states remain non-elemental until an explicit mapping API is introduced.
"""
function state_element(component::Plankton, state::Symbol)
    state in component.states || throw(ArgumentError("Plankton has no state :$state."))
    return state in _ELEMENTAL_STATES ? state : nothing
end

"""Return the conserved-material element represented by a material Pool."""
@inline element(component::Pool) = component.element

"""Return the intrinsic size-structure specification for `component`."""
@inline size_structure(component::Plankton) = component.size_structure
