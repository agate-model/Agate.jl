using Adapt

"""Container for plankton functional-type specifications.

`PFTSpecification` wraps a `NamedTuple` of per-PFT traits and descriptors.
"""
struct PFTSpecification
    data::Any
end

PFTSpecification(; kwargs...) = PFTSpecification((; kwargs...))

@inline pft_has(pft::PFTSpecification, key::Symbol) = hasproperty(pft.data, key)

@inline function pft_get(pft::PFTSpecification, key::Symbol, default=nothing)
    return pft_has(pft, key) ? getproperty(pft.data, key) : default
end

if !hasmethod(Adapt.adapt_structure, Tuple{Any,NamedTuple})
    @inline function Adapt.adapt_structure(to, nt::NamedTuple{names}) where {names}
        return NamedTuple{names}(map(x -> Adapt.adapt(to, x), values(nt)))
    end
end
