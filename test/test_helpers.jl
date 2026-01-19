"""Small test-only helpers.

These exist to keep unit tests independent of specific Oceananigans grid
constructors. In Agate, **grid element type decides precision**, so tests can
use a minimal grid object that exposes `eltype(::grid)` and
`Oceananigans.Architectures.architecture(::grid)`.
"""

import Oceananigans.Architectures: architecture, CPU, GPU

"""A minimal grid stand-in for testing constructor precision/architecture inference.

`Oceananigans` determines the active architecture from the grid. For CPU architectures,
`architecture(grid)` is typically a singleton like `CPU()`, but GPU architectures
carry backend state (e.g. a CUDA backend) and are not nullary-constructible.

This test grid therefore stores an *architecture instance* and returns it directly.
"""
struct DummyGrid{FT,Arch}
    arch::Arch
end

Base.eltype(::DummyGrid{FT,Arch}) where {FT,Arch} = FT
architecture(g::DummyGrid) = g.arch

"""Construct a `DummyGrid` that behaves like an Oceananigans grid."""
dummy_grid(::Type{FT}; arch=CPU()) where {FT<:AbstractFloat} =
    DummyGrid{FT,typeof(arch)}(arch)
