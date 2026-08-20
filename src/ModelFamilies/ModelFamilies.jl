"""Registered model-family identities and their canonical scientific definitions."""
module ModelFamilies

export AbstractModelFamily

"""Registration token for a named biogeochemical model family."""
abstract type AbstractModelFamily end

include("interface.jl")

end
