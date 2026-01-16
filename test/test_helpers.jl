"""Small test-only helpers.

These exist to keep unit tests independent of specific Oceananigans grid
constructors. In Agate, **grid element type decides precision**, so tests can
use a minimal grid object that exposes `eltype(::grid)` and
`Oceananigans.Architectures.architecture(::grid)`.
"""

import Oceananigans.Architectures: architecture, CPU, GPU

"""A minimal grid stand-in for testing constructor precision/architecture inference."""
struct DummyGrid{FT, Arch} end

Base.eltype(::DummyGrid{FT, Arch}) where {FT, Arch} = FT
architecture(::DummyGrid{FT, Arch}) where {FT, Arch} = Arch()

"""Construct a `DummyGrid` that behaves like an Oceananigans grid."""
dummy_grid(::Type{FT}; arch=CPU()) where {FT<:AbstractFloat} = DummyGrid{FT, typeof(arch)}()
