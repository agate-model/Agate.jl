"""Small test-only helpers.

These exist to keep unit tests independent of specific Oceananigans grid
constructors. In Agate, **grid element type decides precision**, so tests can
use a minimal grid object that exposes `eltype(::grid)` and
`Oceananigans.Architectures.architecture(::grid)`.
"""

import OceanBioME
import Oceananigans.Architectures: architecture, CPU

using Agate.Library.Allometry: AllometricParam, ConstantParam, PowerLaw
using Oceananigans.Units: day

const PROCESS_COMPILER_RTOL = 1e-12
process_compiler_isapprox(x, y) =
    isapprox(x, y; rtol=PROCESS_COMPILER_RTOL, atol=10eps(max(abs(x), abs(y))))

"""A minimal grid stand-in for testing constructor precision/architecture inference.

`Oceananigans` determines the active architecture from the grid. For CPU architectures,
`architecture(grid)` is typically a singleton like `CPU()`, but GPU architectures
carry backend state (e.g. a CUDA backend) and are not nullary-constructible.

This test grid therefore stores an *architecture instance* and returns it directly.
"""
struct DummyGrid{T,Arch}
    arch::Arch
end

Base.eltype(::DummyGrid{T,Arch}) where {T,Arch} = T
architecture(g::DummyGrid) = g.arch

"""Construct a `DummyGrid` that behaves like an Oceananigans grid."""
dummy_grid(::Type{T}; arch=CPU()) where {T<:AbstractFloat} = DummyGrid{T,typeof(arch)}(arch)


function authored_nipizd_inputs(::Type{T}=Float32) where {T<:AbstractFloat}
    return (;
        size_structure=(;
            phytoplankton=(diat=T[2, 8],),
            zooplankton=(;
                microzoo=(n=2, min_esd=T(30), max_esd=T(90), splitting=:log_splitting),
            ),
        ),
        scalar_type=T,
        parameters=(;
            maximum_growth_rate=(diat_2=T(1.25 / day),),
            linear_mortality=AllometricParam(
                PowerLaw(); prefactor=T(0.05 / day), exponent=T(-0.1)
            ),
            alpha=ConstantParam(T(0.2 / day)),
        ),
        palatability_matrix=T[0.8 0.2; 0.3 0.7],
        sinking_tracers=(D=T(2.5 / day),),
        open_bottom=false,
    )
end

function nipizd_manifest(
    recipe::Agate.Construction.ProcessModelRecipe;
    grid=OceanBioME.BoxModelGrid(),
    arch=nothing,
    scalar_type=nothing,
)
    _, manifest = Agate.Construction.construct_plus_manifest(
        recipe; grid, arch, scalar_type
    )
    return manifest
end

"""Canonical NiPiZD population realization used by compiler-focused tests."""
default_nipizd_realization() =
    Agate.Models.NiPiZD._population_realization(Agate.Models.NiPiZD.DEFAULT_SIZE_STRUCTURE)

function default_nipizd_layout(::Type{T}=Float64; auxiliary_fields=()) where {T<:Real}
    family = Agate.Models.NiPiZD.NiPiZDFamily()
    realization = default_nipizd_realization()
    return Agate.Configuration.realize_model_layout(
        Agate.ModelFamilies.default_components(family);
        scalar_type=T,
        realization...,
        interaction_roles=(consumers=(:Z,), prey=(:P,)),
        auxiliary_fields,
    )
end
