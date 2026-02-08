"""Factory interfaces.

This module contains the shared base types used by model families to expose
construction-time defaults.

Factories are intentionally minimal: they declare defaults for community structure and
dynamics builders. Group ordering is inferred from the *explicit* ordering of the
`community::NamedTuple` passed to `Construction.construct_factory`.
"""
module Factories

export AbstractBGCFactory

"""Abstract supertype for biogeochemical model factories."""
abstract type AbstractBGCFactory end

end # module
