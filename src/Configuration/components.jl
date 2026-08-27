"""Population state described by intrinsic component properties.

Ecological function is supplied by process participation. `states` maps each
prognostic state identity to the conserved-material currency it represents,
while `size_structure` describes the shared ecological class realization.
A one-state population is authored with its currency as the positional argument.
"""
struct Population{ST<:NamedTuple,S}
    states::ST
    size_structure::S

    function Population(states::ST, size_structure::S) where {ST<:NamedTuple,S}
        isempty(states) && throw(ArgumentError("Population must define at least one state."))
        all(value -> value isa Symbol, values(states)) || throw(
            ArgumentError("Population state currencies must be Symbols."),
        )
        return new{ST,S}(states, size_structure)
    end
end

function Population(currency; size_structure=nothing)
    currency isa Symbol || throw(ArgumentError("Population currency must be a Symbol."))
    return Population(NamedTuple{(currency,)}((currency,)), size_structure)
end

function Population(; states, size_structure=nothing)
    states isa NamedTuple || throw(
        ArgumentError("Population states must be a NamedTuple mapping state names to currencies."),
    )
    return Population(states, size_structure)
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

"""Construct a setup-time reference to one population prognostic state."""
@inline population_state(population::Symbol, state::Symbol) = PopulationStateRef(population, state)

"""Return the prognostic state-to-currency mapping for a population."""
@inline states(component::Population) = component.states

"""Return the currency represented by one named population state."""
function state_currency(component::Population, state::Symbol)
    hasproperty(component.states, state) || throw(
        ArgumentError("Population has no state :$state."),
    )
    return getproperty(component.states, state)
end

"""Return the conserved-material currency represented by a single-state component."""
@inline currency(component::Pool) = component.currency
function currency(component::Population)
    length(component.states) == 1 || throw(
        ArgumentError(
            "multi-state Population has no unique currency; use state_currency(population, state)",
        ),
    )
    return only(values(component.states))
end

"""Return the intrinsic size-structure specification for `component`."""
@inline size_structure(component::Union{Population,Pool}) = component.size_structure

