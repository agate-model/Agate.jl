using Adapt

if !hasmethod(Adapt.adapt_structure, Tuple{Any,NamedTuple})
    @inline function Adapt.adapt_structure(to, nt::NamedTuple{Names}) where {Names}
        return NamedTuple{Names}(map(x -> Adapt.adapt(to, x), values(nt)))
    end
end
