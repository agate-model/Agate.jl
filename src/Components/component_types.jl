"""Plankton state described by intrinsic component properties.

Ecological function is supplied by process participation. `states` declares each prognostic
state together with the conserved Element it represents, using `nothing` for non-elemental
states. `reference_state` identifies the state used as the biological reference basis.

For example, `states=(biomass=:carbon, chlorophyll=nothing)` declares a carbon inventory named
`:biomass` plus a non-elemental chlorophyll state. Element identity is therefore independent of
state spelling and is queried centrally through [`state_element`](@ref). Each conserved Element may
currently be represented by at most one state. `composition` maps implicitly represented Elements
to model parameter names giving their amount per unit reference-state Element; explicit states
continue to represent independent, variable inventories.
"""
struct Plankton{StateElements<:NamedTuple,SizeStructure,Composition<:NamedTuple}
    state_elements::StateElements
    reference_state::Symbol
    size_structure::SizeStructure
    composition::Composition

    function Plankton(
        state_elements::StateElements, reference_state::Symbol, size_structure::SizeStructure,
        composition::Composition,
    ) where {StateElements<:NamedTuple,SizeStructure,Composition<:NamedTuple}
        isempty(state_elements) && throw(ArgumentError("Plankton must define at least one state."))
        all(value -> isnothing(value) || value isa Symbol, values(state_elements)) || throw(
            ArgumentError("Plankton state Elements must be Symbols or `nothing`."),
        )
        elemental = Tuple(value for value in values(state_elements) if !isnothing(value))
        length(unique(elemental)) == length(elemental) || throw(
            ArgumentError("Plankton states must represent distinct Elements."),
        )
        reference_state in keys(state_elements) || throw(
            ArgumentError("Plankton `reference_state` must be one of its declared states."),
        )
        reference_element = getproperty(state_elements, reference_state)
        all(value -> value isa Symbol, values(composition)) || throw(
            ArgumentError("Plankton `composition` must map Elements to parameter Symbols."),
        )
        isempty(composition) || reference_element isa Symbol || throw(
            ArgumentError("Plankton with implicit composition requires an elemental reference state."),
        )
        for element in keys(composition)
            element === reference_element && throw(ArgumentError(
                "Plankton `composition` cannot repeat its reference element :$element.",
            ))
            element in elemental && throw(ArgumentError(
                "Plankton `composition` element :$element already has a prognostic state.",
            ))
        end
        return new{StateElements,SizeStructure,Composition}(
            state_elements, reference_state, size_structure, composition,
        )
    end
end

function Plankton(; states, reference_state, size_structure=nothing, composition=(;))
    states isa NamedTuple || throw(ArgumentError(
        "Plankton `states` must be a NamedTuple mapping state names to Elements or `nothing`.",
    ))
    reference_state isa Symbol || throw(
        ArgumentError("Plankton `reference_state` must be a Symbol."),
    )
    composition isa NamedTuple || throw(ArgumentError(
        "Plankton `composition` must be a NamedTuple mapping Elements to parameter Symbols.",
    ))
    names = Tuple(sort!(collect(keys(composition)); by=String))
    canonical = NamedTuple{names}(Tuple(getproperty(composition, name) for name in names))
    return Plankton(states, reference_state, size_structure, canonical)
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
@inline states(component::Plankton) = keys(component.state_elements)

"""Return the plankton state used as its biological reference basis."""
@inline reference_state(component::Plankton) = component.reference_state

"""Return all prognostic plankton states other than the reference state."""
@inline variable_states(component::Plankton) = Tuple(
    state for state in states(component) if state !== component.reference_state
)

"""Return the conserved Element represented by `state`, or `nothing` for a non-elemental state.

Element identity is the additive elemental inventory represented by the state's numerical value,
not the chemical composition of the represented material. State identity and Element identity are
independent: any declared state name may map explicitly to an Element Symbol.
"""
function state_element(component::Plankton, state::Symbol)
    hasproperty(component.state_elements, state) || throw(
        ArgumentError("Plankton has no state :$state."),
    )
    return getproperty(component.state_elements, state)
end

"""Return the conserved-material element represented by a material Pool."""
@inline element(component::Pool) = component.element

"""Return the intrinsic size-structure specification for `component`."""
@inline size_structure(component::Plankton) = component.size_structure

"""Implicit elemental ratios (amount per unit reference element), keyed by Element."""
@inline composition(component::Plankton) = component.composition
