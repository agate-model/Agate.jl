# paper/GPU/double_gyre_gpu_online_NiPiZD_noCC_noFW.jl
#
# Double gyre (GPU, online physics) + Agate NiPiZD (Z1, Z2, P1, P2, D, N) with light attenuation.
# Seasonal forcing kept; NO warming; NO freshwater flux; max Δt = 180 minutes.
#
# Run:
#   julia --project=paper/GPU paper/GPU/double_gyre_gpu_online_NiPiZD_noCC_noFW.jl 70 30 15
#
# Args:
#   1: stop_years (default 70)
#   2: Δt0_minutes (default 30)
#   3: out_interval_days (default 15; set 0 to disable output)

using Oceananigans
using Oceananigans.Units: minute, hour, day
using Oceananigans.Advection: UpwindBiased
using Oceananigans.OutputWriters: JLD2Writer
using Printf

using CUDA
CUDA.allowscalar(false)
using JLD2

using Agate
using Agate.Library.Light
using Agate.Models: NiPiZD
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers
using Adapt

if !CUDA.functional()
    error("CUDA is not functional. Check nvidia-smi and CUDA.jl setup inside the container.")
end

const FT   = Float32
const year = 365day
const arch = GPU()

# -------------------------
# Domain + grid (1°-like)
# -------------------------
const Lx = 3180e3
const Ly = 3180e3
const H  = 4000.0
const Nx = 30
const Ny = 30
const Nz = 30

function couespel_z_faces(::Type{FT}; H=4000.0, Nz=30) where {FT}
    Δz1 = collect(range(FT(10.0),  FT(14.0);  length=10))
    Δz2 = collect(range(FT(20.0),  FT(56.0);  length=10))
    Δz3 = collect(range(FT(100.0), FT(500.0); length=10))

    Δz1 .*= FT(120.0) / sum(Δz1)
    Δz2 .*= FT(380.0) / sum(Δz2)
    Δz3 .*= FT(H - 500.0) / sum(Δz3)

    z_faces = zeros(FT, Nz+1)
    z_faces[1] = -FT(H)
    for k in 1:Nz
        z_faces[k+1] = z_faces[k] + (k <= 10 ? Δz1[k] : k <= 20 ? Δz2[k-10] : Δz3[k-20])
    end
    z_faces[end] = FT(0)
    return z_faces
end

z_faces = couespel_z_faces(FT; H=H, Nz=Nz)

grid = RectilinearGrid(arch, FT;
                       size=(Nx, Ny, Nz),
                       x = (-FT(Lx/2), FT(Lx/2)),
                       y = (-FT(Ly/2), FT(Ly/2)),
                       z = z_faces,
                       topology = (Bounded, Bounded, Bounded))

# -------------------------
# Coriolis + buoyancy
# -------------------------
coriolis = BetaPlane(; latitude=FT(30), radius=FT(6.371e6), rotation_rate=FT(7.292115e-5))

const α  = FT(2.0e-4)
const βS = FT(7.6e-4)
eos = LinearEquationOfState(; thermal_expansion=α, haline_contraction=βS)
buoyancy = SeawaterBuoyancy(; equation_of_state=eos)

# -------------------------
# Mixing / closure
# -------------------------
const νh    = FT(1.0e5)
const κh    = FT(1000.0)
const νz_bg = FT(1.0e-4)
const κz_bg = FT(1.0e-5)

function build_closure()
    Az = HorizontalScalarDiffusivity(; ν=νh, κ=κh)
    Kz = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(); ν=νz_bg, κ=κz_bg)
    return (Az, Kz)
end

# -------------------------
# Seasonal forcing only (NO warming, NO freshwater flux)
# -------------------------
const ρ0 = FT(1026.0)
const τ0       = FT(0.10)    # N/m²
const τ_season = FT(0.015)
const h_surface  = FT(50.0)  # m
const τT_restore = 30day
const S0 = FT(35.0)

@inline ηy(y) = (y + FT(Ly/2)) / FT(Ly)

@inline function τx_surface(y, t)
    seasonal = τ0 + τ_season * FT(sin(2π * t / year))
    return -seasonal * FT(cos(2π * ηy(y)))
end

@inline function T_air(y, t)
    η = ηy(y)
    Tmean   = FT(15.0) + FT(10.0) * FT(cos(π * η))
    Tseason = FT(5.0)  * FT(sin(2π * t / year))
    return Tmean + Tseason
end

@inline surface_mask(x, y, z) = z > -h_surface ? FT(1) : FT(0)

# -------------------------
# Light forcing + attenuation (Agate)
# -------------------------
@inline function PAR0(x, y, t)
    seasonal = 1 - cos((t + 15day) * 2π / year)
    gaussian = exp(-((mod(t, year) - 200day) / (50day))^2)
    pulse    = 1 / (1 + 0.2 * gaussian) + 2
    return FT(60) * FT(seasonal) * FT(pulse)
end

@inline PAR_f(x, y, z, t) = PAR0(x, y, t) * exp(FT(0.2) * z)

light_attenuation = FunctionFieldPAR(; grid, PAR_f=PAR_f)

# -------------------------
# Build NiPiZD biogeochem (GPU) — robust to construct() returning instance vs builder
# -------------------------
bgc_obj = NiPiZD.construct()

bgc_cpu = try
    bgc_obj()
catch
    bgc_obj
end

bgc_tracers_any = required_biogeochemical_tracers(bgc_cpu)
bgc_tracers = bgc_tracers_any isa Symbol ? (bgc_tracers_any,) : Tuple(bgc_tracers_any)

bgc_gpu = Adapt.adapt(CuArray, bgc_cpu)
biogeochemistry = Biogeochemistry(bgc_gpu; light_attenuation)

# -------------------------
# Boundary conditions (closed basin)
# -------------------------
u_bcs = FieldBoundaryConditions(
    top    = FluxBoundaryCondition((x, y, t) -> τx_surface(y, t) / ρ0),
    bottom = FluxBoundaryCondition(FT(0)),
    south  = FluxBoundaryCondition(FT(0)),
    north  = FluxBoundaryCondition(FT(0))
)

v_bcs = FieldBoundaryConditions(
    top    = FluxBoundaryCondition(FT(0)),
    bottom = FluxBoundaryCondition(FT(0)),
    west   = FluxBoundaryCondition(FT(0)),
    east   = FluxBoundaryCondition(FT(0))
)

T_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(FT(0)),
                                bottom=FluxBoundaryCondition(FT(0)),
                                west=FluxBoundaryCondition(FT(0)),
                                east=FluxBoundaryCondition(FT(0)),
                                south=FluxBoundaryCondition(FT(0)),
                                north=FluxBoundaryCondition(FT(0)))

S_bcs = T_bcs

# -------------------------
# Model builder
# -------------------------
function build_model()
    closure = build_closure()
    adv = UpwindBiased(order=3)

    T_restoring = Relaxation(; rate   = 1 / τT_restore,
                               mask   = surface_mask,
                               target = (x, y, z, t) -> T_air(y, t))

    all_tracers = Tuple(unique([:T, :S, bgc_tracers...]))

    model = HydrostaticFreeSurfaceModel(
        grid = grid,
        free_surface = ImplicitFreeSurface(),
        coriolis = coriolis,
        buoyancy = buoyancy,
        tracers = all_tracers,
        closure = closure,
        boundary_conditions = (u=u_bcs, v=v_bcs, T=T_bcs, S=S_bcs),

        momentum_advection = adv,
        tracer_advection   = adv,

        timestepper = :SplitRungeKutta3,

        forcing = (T = T_restoring,),
        biogeochemistry = biogeochemistry
    )

    return model
end

# -------------------------
# Initial conditions
# -------------------------
@inline T_init(z) = FT(20.0) + (FT(5.0) - FT(20.0)) * (-z / FT(H))

function set_initial_conditions!(model)
    # physics tracers
    set!(model,
         u = FT(0), v = FT(0),
         T = (x, y, z) -> T_init(z),
         S = S0)

    # biogeochem tracers: set fields directly (no dynamic keywords)
    defaults = Dict{Symbol, FT}(
        :P1 => FT(0.01), :P2 => FT(0.01),
        :Z1 => FT(0.05), :Z2 => FT(0.05),
        :N  => FT(7.0),
        :D  => FT(1.0)
    )

    for (name, val) in defaults
        if haskey(model.tracers, name)
            set!(model.tracers[name], val)
        end
    end

    return nothing
end

# -------------------------
# Run
# -------------------------
function run_simulation(; stop_years::Float64=70.0,
                          Δt0 = 30minute,
                          out_interval_days::Int=15,
                          filename::String="double_gyre_NiPiZD_noCC_noFW")

    model = build_model()
    set_initial_conditions!(model)

    simulation = Simulation(model; Δt=Δt0, stop_time=stop_years * year)

    wizard = TimeStepWizard(cfl=0.6, max_change=1.2, min_change=0.5, max_Δt=180minute)
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

    function progress(sim)
        c = sim.model.clock
        @printf("iter: %d | time: %.4f yr | Δt: %.2f min\n",
                c.iteration, c.time / year, sim.Δt / minute)
        return nothing
    end
    simulation.callbacks[:progress] = Callback(progress, IterationInterval(200))

    if out_interval_days > 0
        simulation.output_writers[:profiles] = JLD2Writer(
            model,
            model.tracers;
            filename = "$(filename).jld2",
            schedule = TimeInterval(out_interval_days * day),
            overwrite_existing = true
        )
        @info "Writing JLD2 outputs to $(filename).jld2 every $(out_interval_days) days"
    else
        @info "Output disabled (out_interval_days <= 0)."
    end

    run!(simulation)
    return nothing
end

# -------------------------
# CLI
# -------------------------
stop_years        = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 70.0
Δt0_minutes       = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 30.0
out_interval_days = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 15

run_simulation(; stop_years=stop_years,
               Δt0=Δt0_minutes * minute,
               out_interval_days=out_interval_days)