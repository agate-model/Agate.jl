"""Small test-only helpers.

These exist to keep unit tests independent of specific Oceananigans grid
constructors. In Agate, **grid element type decides precision**, so tests can
use a minimal grid object that exposes `eltype(::grid)` and
`Oceananigans.Architectures.architecture(::grid)`.
"""

import Oceananigans.Architectures: architecture, CPU

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

const MULTI_NUTRIENT_COUPLINGS = (
    Agate.Tendencies.nutrient_coupling(
        :DIN,
        :half_saturation_DIN;
        stoichiometry=:nitrogen_to_carbon,
        remineralization=((:DON, :organic_remineralization), (:PON, :organic_remineralization)),
    ),
    Agate.Tendencies.nutrient_coupling(
        :PO4,
        :half_saturation_PO4;
        stoichiometry=:phosphorus_to_carbon,
        remineralization=((:DOP, :organic_remineralization), (:POP, :organic_remineralization)),
    ),
)

function multi_nutrient_config(limitation)
    return Agate.Tendencies.TendencyConfig(;
        growth=:smith,
        organic_cycling=:dom_pom,
        nutrient_limitation=limitation,
        nutrients=MULTI_NUTRIENT_COUPLINGS,
    )
end

const MULTI_NUTRIENT_LIEBIG = multi_nutrient_config(:liebig)
const MULTI_NUTRIENT_FRANK = multi_nutrient_config(Agate.Library.Nutrients.FrankTNorm(50))

function multi_nutrient_test_model()
    interactions = (consumer_global=Int[], prey_global=Int[], global_to_prey=[0])
    parameters = (
        organic_remineralization=0.1213 / 86400,
        nitrogen_to_carbon=16 / 106,
        phosphorus_to_carbon=1 / 106,
        DOM_POM_fractionation=0.5,
        linear_mortality=[8e-7],
        quadratic_mortality=[0.0],
        maximum_growth_rate=[2 / 86400],
        half_saturation_DIN=[0.5],
        half_saturation_PO4=[0.5],
        alpha=[0.1 / 86400],
        maximum_predation_rate=[0.0],
        interactions=interactions,
    )
    config = MULTI_NUTRIENT_LIEBIG
    organic(target, fraction, stoichiometry=:one) = Agate.Tendencies.organic_matter_tendency(
        config; target, remineralization=:organic_remineralization, fraction, stoichiometry
    )
    tracers = (
        DIC=Agate.Tendencies.inorganic_tendency(
            config;
            target=:DIC,
            remineralization=((:DOC, :organic_remineralization), (:POC, :organic_remineralization)),
            stoichiometry=:one,
        ),
        DIN=Agate.Tendencies.inorganic_tendency(config; target=:DIN),
        PO4=Agate.Tendencies.inorganic_tendency(config; target=:PO4),
        DOC=organic(:DOC, :DOM),
        POC=organic(:POC, :POM),
        DON=organic(:DON, :DOM, :nitrogen_to_carbon),
        PON=organic(:PON, :POM, :nitrogen_to_carbon),
        DOP=organic(:DOP, :DOM, :phosphorus_to_carbon),
        POP=organic(:POP, :POM, :phosphorus_to_carbon),
        P_1=Agate.Tendencies.phytoplankton_tendency(config; plankton_idx=1),
    )
    community = (group_symbols=(:P,), group_indices=(P=1:1,))
    tracer_index = Agate.Runtime.build_tracer_index(
        community, keys(tracers), (:PAR,); n_biogeochem_tracers=9
    )
    factory = Agate.Construction.define_tracer_functions(
        parameters, tracers; tracer_index
    )
    return factory(parameters)
end
