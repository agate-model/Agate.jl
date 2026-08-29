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

# Shared quota-model fixture used by authoring, compilation, and AD tests.
quota_components() = (
    DIC=Agate.Configuration.Pool(:carbon),
    DIN=Agate.Configuration.Pool(:nitrogen),
    PO4=Agate.Configuration.Pool(:phosphorus),
    P=Agate.Configuration.Plankton(;
        states=(:carbon, :nitrogen, :phosphorus),
        reference_state=:carbon,
        size_structure=[1.0, 2.0],
    ),
)

quota_response(state, minimum, maximum) = Agate.Processes.QuotaResponse(
    Agate.Processes.NormalizedDroop();
    variable_state=state,
    bindings=(minimum_quota=minimum, maximum_quota=maximum),
)

quota_uptake(state, resource, bindings) = Agate.Processes.NutrientUptake(
    Agate.Processes.QuotaRegulatedMonod();
    plankton=:P,
    target_state=state,
    resource=resource,
    bindings=bindings,
)

function quota_processes()
    responses = (
        nitrogen=quota_response(
            :nitrogen, :minimum_nitrogen_quota, :maximum_nitrogen_quota
        ),
        phosphorus=quota_response(
            :phosphorus, :minimum_phosphorus_quota, :maximum_phosphorus_quota
        ),
    )
    growth = Agate.Processes.Growth(;
        plankton=:P,
        reference_resource=:DIC,
        bindings=(maximum_rate=:maximum_growth_rate,),
        factors=(
            light=Agate.Processes.Light(
                Agate.Processes.Smith(); driver=:PAR,
                bindings=(alpha=:photosynthetic_slope,),
            ),
            nutrients=Agate.Processes.NutrientLimitation(
                Agate.Processes.Liebig(); responses=responses
            ),
        ),
    )
    nitrogen_uptake = quota_uptake(:nitrogen, :DIN, (
        maximum_rate=:maximum_nitrogen_uptake,
        half_saturation=:nitrogen_half_saturation,
        minimum_quota=:minimum_nitrogen_quota,
        maximum_quota=:maximum_nitrogen_quota,
        hill=:nitrogen_uptake_hill,
    ))
    phosphorus_uptake = quota_uptake(:phosphorus, :PO4, (
        maximum_rate=:maximum_phosphorus_uptake,
        half_saturation=:phosphorus_half_saturation,
        minimum_quota=:minimum_phosphorus_quota,
        maximum_quota=:maximum_phosphorus_quota,
        hill=:phosphorus_uptake_hill,
    ))
    return (; growth, nitrogen_uptake, phosphorus_uptake)
end

quota_parameters() = (
    maximum_growth_rate=Agate.Parameters.Parameter(0.5),
    photosynthetic_slope=Agate.Parameters.Parameter(0.05),
    minimum_nitrogen_quota=Agate.Parameters.Parameter(0.05),
    maximum_nitrogen_quota=Agate.Parameters.Parameter(0.2),
    minimum_phosphorus_quota=Agate.Parameters.Parameter(0.005),
    maximum_phosphorus_quota=Agate.Parameters.Parameter(0.02),
    maximum_nitrogen_uptake=Agate.Parameters.Parameter(0.1),
    nitrogen_half_saturation=Agate.Parameters.Parameter(0.2),
    nitrogen_uptake_hill=Agate.Parameters.Parameter(2.0),
    maximum_phosphorus_uptake=Agate.Parameters.Parameter(0.01),
    phosphorus_half_saturation=Agate.Parameters.Parameter(0.02),
    phosphorus_uptake_hill=Agate.Parameters.Parameter(2.0),
)

quota_definition() = Agate.Processes.ModelDefinition(;
    components=quota_components(),
    processes=quota_processes(),
    parameters=quota_parameters(),
)

nipizd_named_size_structure() = (;
    phytoplankton=(diat=[2.0, 5.0, 10.0], dino=[8.0, 20.0]),
    zooplankton=(microzoo=[30.0, 60.0], mesozoo=[100.0]),
)

nipizd_test_state(::Type{T}=Float64) where {T<:AbstractFloat} = (;
    N=T(7), D=T(1), P_1=T(0.01), P_2=T(0.01), Z_1=T(0.05), Z_2=T(0.05),
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

function test_tendencies(model, state::NamedTuple, expected::NamedTuple)
    names = Agate.Introspection.tracer_names(model)
    args = (0.0, 0.0, 0.0, 0.0, Tuple(getproperty(state, name) for name in names)...)
    for (tracer, tendency) in pairs(expected)
        @test model(Val(tracer), args...) ≈ tendency
    end
end

function model_tendencies(bgc, args; tracers=Tuple(Agate.Introspection.tracer_names(bgc)))
    return NamedTuple{tracers}(Tuple(bgc(Val(tracer), args...) for tracer in tracers))
end

function authored_nipizd_inputs(::Type{T}=Float32) where {T<:AbstractFloat}
    return (;
        size_structure=(;
            phytoplankton=(diat=T[2, 8],),
            zooplankton=(;
                microzoo=(n=2, min_esd=T(30), max_esd=T(90), spacing=:log),
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
    recipe::Agate.Construction.ModelRecipe;
    grid=OceanBioME.BoxModelGrid(),
    arch=nothing,
    scalar_type=nothing,
)
    _, manifest = Agate.Construction.construct_plus_manifest(
        recipe; grid, arch, scalar_type
    )
    return manifest
end
