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

function nipizd_recipe_manifest(; kwargs...)
    inputs = Agate.Models.NiPiZD._construction_inputs(; kwargs...)
    recipe = Agate.Construction.capture_model_recipe(
        inputs.factory; inputs.recipe_kwargs...
    )
    _, manifest = Agate.Construction.construct_factory_plus_manifest(
        inputs.factory; inputs.kwargs...
    )
    return recipe, manifest
end

function nipizd_manifest(
    recipe::Agate.Construction.ModelRecipe;
    grid=OceanBioME.BoxModelGrid(),
    arch=nothing,
    scalar_type=nothing,
)
    inputs = Agate.Models.NiPiZD._recipe_construction_inputs(
        recipe; grid, arch, scalar_type
    )
    _, manifest = Agate.Construction.construct_factory_plus_manifest(
        inputs.factory; inputs.kwargs...
    )
    return manifest
end

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
