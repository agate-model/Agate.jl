using Oceananigans
using Oceananigans.Units: minute, hour, day
using Oceananigans.Advection: UpwindBiased
using Oceananigans.OutputWriters: JLD2Writer
using Printf

using CUDA
CUDA.allowscalar(false)

using JLD2

if !CUDA.functional()
    error("CUDA is not functional in this Julia session. Check GPU passthrough (nvidia-smi) and CUDA.jl setup.")
end

const FT = Float32
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
# Seasonal forcing 
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
# Boundary conditions (closed basin, free-slip-ish)
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
# Model
# -------------------------
function build_model()
    free_surface = ImplicitFreeSurface()
    closure = build_closure()

    T_restoring = Relaxation(; rate   = 1 / τT_restore,
                               mask   = surface_mask,
                               target = (x, y, z, t) -> T_air(y, t))

    adv = UpwindBiased(order=3)

    return HydrostaticFreeSurfaceModel(
        grid = grid,
        free_surface = free_surface,
        coriolis = coriolis,
        buoyancy = buoyancy,
        tracers = (:T, :S),
        closure = closure,
        boundary_conditions = (u=u_bcs, v=v_bcs, T=T_bcs, S=S_bcs),

        momentum_advection = adv,
        tracer_advection   = adv,

        timestepper = :SplitRungeKutta3,

        forcing = (T = T_restoring,)  
    )
end

# -------------------------
# Initial conditions
# -------------------------
@inline T_init(z) = FT(20.0) + (FT(5.0) - FT(20.0)) * (-z / FT(H))
@inline S_init(z) = S0

function set_initial_conditions!(model)
    set!(model,
         u = FT(0), v = FT(0),
         T = (x, y, z) -> T_init(z),
         S = (x, y, z) -> S_init(z))
end

# -------------------------
# Run
# -------------------------
function run_simulation(; stop_years::Float64=70.0,
                          Δt0 = 30minute,
                          filename::String="double_gyre_online_noCC_noFW",
                          out_interval_days::Int=180)

    model = build_model()
    set_initial_conditions!(model)

    simulation = Simulation(model; Δt=Δt0, stop_time=stop_years * year)

    # Allow Δt up to 180 minutes (3 hours)
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
        # Like your working pattern:
        # - pass a NamedTuple of fields
        # - write at TimeInterval(out_interval_days * day)
        fields_to_write = merge(model.tracers, (; η = model.free_surface.η))

        simulation.output_writers[:profiles] = JLD2Writer(
            model,
            fields_to_write;
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
# Usage:
#   julia --project=paper/GPU paper/GPU/double_gyre_gpu_online_noCC_noFW.jl
#   julia --project=paper/GPU paper/GPU/double_gyre_gpu_online_noCC_noFW.jl 70 30 180
#
# Args:
#   1: stop_years (default 70)
#   2: Δt0_minutes (default 30)
#   3: out_interval_days (default 180; set 0 to disable output)
stop_years = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 70.0
Δt0_minutes = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 30.0
out_interval_days = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 15

run_simulation(; stop_years=stop_years,
               Δt0=Δt0_minutes * minute,
               out_interval_days=out_interval_days)