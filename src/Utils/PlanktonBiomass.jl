
struct PlanktonBiomass{FT <: Float64} <: AbstractVector{FT}
    P1::FT
    P2::FT
    Z1::FT
    Z2::FT
end

# Array interface
Base.size(::PlanktonBiomass) = (4,)
Base.getindex(b::PlanktonBiomass, i::Int) = getfield(b, i)
Base.setindex!(b::PlanktonBiomass, v, i::Int) = setfield!(b, i, v)

# Broadcasting
Base.BroadcastStyle(::Type{<:PlanktonBiomass}) = Broadcast.DefaultArrayStyle{1}()
# Base.BroadcastStyle(::Type{<:PlanktonBiomass}) = Broadcast.ArrayStyle{PlanktonBiomass}()
Base.similar(b::PlanktonBiomass{FT}) where {FT} = PlanktonBiomass(zero(FT), zero(FT), zero(FT), zero(FT))

# Adjoint for linear algebra
using LinearAlgebra
Base.adjoint(pb::PlanktonBiomass) = Adjoint(pb) 

# GPU compatibility
Adapt.@adapt_structure PlanktonBiomass