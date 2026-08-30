using Adapt

if !hasmethod(Adapt.adapt_structure, Tuple{Any,NamedTuple})
    @inline function Adapt.adapt_structure(to, nt::NamedTuple{names}) where {names}
        return NamedTuple{names}(map(x -> Adapt.adapt(to, x), values(nt)))
    end
end
