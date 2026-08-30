"""Interface for named model-family identities and scientific definition hooks."""
module ModelFamilies

export AbstractModelFamily

"""Identity token for a named biogeochemical model family."""
abstract type AbstractModelFamily end

include("interface.jl")

end
