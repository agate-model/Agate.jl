"""Small test-only helpers shared across the behavioral suite."""

import OceanBioME
import Oceananigans.Architectures: architecture, CPU

using Agate.Library.Allometry: AllometricParam, ConstantParam, PowerLaw
using Oceananigans.Units: day

const PROCESS_COMPILER_RTOL = 1e-12
process_compiler_isapprox(x, y) =
    isapprox(x, y; rtol=PROCESS_COMPILER_RTOL, atol=10eps(max(abs(x), abs(y))))

struct DummyGrid{T,Arch}
    arch::Arch
end

Base.eltype(::DummyGrid{T,Arch}) where {T,Arch} = T
architecture(g::DummyGrid) = g.arch

dummy_grid(::Type{T}; arch=CPU()) where {T<:AbstractFloat} = DummyGrid{T,typeof(arch)}(arch)

function argument_error_message(f)
    error = try
        f()
        nothing
    catch caught
        caught
    end
    @test error isa ArgumentError
    return error isa Exception ? sprint(showerror, error) : ""
end

canonicalization_error_message(definition) =
    argument_error_message(() -> Agate.Processes.canonicalize_model(definition))

nipizd_named_size_structure() = (;
    phytoplankton=(diat=[2.0, 5.0, 10.0], dino=[8.0, 20.0]),
    zooplankton=(microzoo=[30.0, 60.0], mesozoo=[100.0]),
)

nipizd_test_state(::Type{T}=Float64) where {T<:AbstractFloat} = (;
    N=T(7), D=T(1), Z_1=T(0.05), Z_2=T(0.05), P_1=T(0.01), P_2=T(0.01),
)

nipizd_u0(::Type{T}=Float64) where {T<:AbstractFloat} = collect(values(nipizd_test_state(T)))

function nipizd_runtime_args(::Type{T}=Float64; PAR=100) where {T<:AbstractFloat}
    return (0, 0, 0, 0, values(nipizd_test_state(T))..., T(PAR))
end

function nipizd_growth_fixture(; mu=0.7 / day, kwargs...)
    return Agate.Models.NiPiZD.construct(;
        parameters=(; maximum_growth_rate=(P_1=mu, P_2=mu)), kwargs...
    )
end

function model_tendencies(bgc, args; tracers=Tuple(Agate.Introspection.tracer_names(bgc)))
    return NamedTuple{tracers}(Tuple(bgc(Val(tracer), args...) for tracer in tracers))
end

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
