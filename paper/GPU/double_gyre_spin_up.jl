# double_gyre_physics_spinup_offline_year.jl
#
# Physics-only double gyre spin-up:
#   - years 1:69  : physics only, no archive
#   - year 70     : physics only, archive u, v, w, T, S for offline cyclic BGC forcing
#
# Outputs:
#   1) <prefix>_cycle_init.jld2     : state at start of archived year (recommended BGC initial state)
#   2) <prefix>_offline_year.jld2   : archived final year for cyclic offline forcing
#   3) final checkpoint             : optional terminal restart for further physics-only continuation
#
# Run:
#   julia --project=. double_gyre_physics_spinup_offline_year.jl 70 45 5
#
# Args:
#   1: total_years        (default 70)
#   2: Δt0_minutes        (default 45)
#   3: archive_days       (default 5)

using Oceananigans
using Oceananigans.Units: minute, hour, day
using Oceananigans.Advection: UpwindBiased
using Oceananigans.OutputWriters: JLD2Writer, Checkpointer
using Oceananigans.Grids: Center, node
using Printf
using CUDA
CUDA.allowscalar(false)
using JLD2

if !CUDA.functional()
    error("CUDA is not functional. Check nvidia-smi and CUDA.jl setup.")
end

const FT   = Float32
const year = 365day
const arch = GPU()

# -------------------------
# Domain + grid
# -------------------------
const Lx = 3180e3
const Ly = 3180e3
const H  = 4000.0
const Nx = 60
const Ny = 60
const Nz = 30

function couespel_z_faces(::Type{FT}; H=4000.0, Nz=30) where {FT}
    Δz1 = collect(range(FT(10.0),  FT(14.0);  length=10))
    Δz2 = collect(range(FT(20.0),  FT(56.0);  length=10))
    Δz3 = collect(range(FT(100.0), FT(500.0); length=10))

    Δz1 .*= FT(120.0) / sum(Δz1)
    Δz2 .*= FT(380.0) / sum(Δz2)
    Δz3 .*= FT(H - 500.0) / sum(Δz3)

    Δz = vcat(Δz3, Δz2, Δz1)

    z_faces = zeros(FT, Nz + 1)
    z_faces[1] = -FT(H)
    for k in 1:Nz
        z_faces[k + 1] = z_faces[k] + Δz[k]
    end
    z_faces[end] = FT(0)
    return z_faces
end

z_faces = couespel_z_faces(FT; H=H, Nz=Nz)

grid = RectilinearGrid(arch, FT;
                       size=(Nx, Ny, Nz),
                       x = (-FT(Lx / 2), FT(Lx / 2)),
                       y = (-FT(Ly / 2), FT(Ly / 2)),
                       z = z_faces,
                       topology = (Bounded, Bounded, Bounded))

# -------------------------
# Coriolis + buoyancy
# -------------------------
coriolis = BetaPlane(; latitude=FT(30),
                     radius=FT(6.371e6),
                     rotation_rate=FT(7.292115e-5))

const α  = FT(2.0e-4)
const βS = FT(7.6e-4)

eos = LinearEquationOfState(; thermal_expansion=α,
                              haline_contraction=βS)

buoyancy = SeawaterBuoyancy(; equation_of_state=eos)

# -------------------------
# Mixing / closure
# -------------------------
const νh    = FT(1.0e5)
const κh    = FT(3.0e3)
const νz_bg = FT(1.0e-4)
const κz_bg = FT(1.0e-5)

function build_closure()
    Az = HorizontalScalarDiffusivity(; ν=νh, κ=κh)
    Kz = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization();
                                   ν=νz_bg, κ=κz_bg)
    return (Az, Kz)
end

# -------------------------
# Seasonal forcing only
# -------------------------
const ρ0         = FT(1026.0)
const τ0         = FT(0.10)
const τ_season   = FT(0.015)
const h_surface  = FT(50.0)
const τT_restore = 30day
const S0         = FT(35.0)

@inline ηy(y) = (y + FT(Ly / 2)) / FT(Ly)

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

# -------------------------
# Boundary conditions
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
@inline function T_restoring_discrete(i, j, k, grid, clock, fields)
    x, y, z = node(i, j, k, grid, Center(), Center(), Center())
    T = @inbounds fields.T[i, j, k]
    return ifelse(z > -h_surface,
                  (T_air(y, clock.time) - T) / τT_restore,
                  FT(0))
end

function build_model()
    closure      = build_closure()
    momentum_adv = UpwindBiased(order=3)
    tracer_adv   = UpwindBiased(order=1)

    T_restoring = Forcing(T_restoring_discrete; discrete_form=true)

    model = HydrostaticFreeSurfaceModel(
        grid = grid,
        free_surface = ImplicitFreeSurface(),
        coriolis = coriolis,
        buoyancy = buoyancy,
        tracers = (:T, :S),
        closure = closure,
        boundary_conditions = (u=u_bcs, v=v_bcs, T=T_bcs, S=S_bcs),
        momentum_advection = momentum_adv,
        tracer_advection   = tracer_adv,
        timestepper = :SplitRungeKutta3,
        forcing = (T = T_restoring,)
    )

    return model
end

# -------------------------
# Initial conditions
# -------------------------
@inline T_init(z) = FT(20.0) + (FT(5.0) - FT(20.0)) * (-z / FT(H))

function set_initial_conditions!(model)
    set!(model,
         u = FT(0),
         v = FT(0),
         T = (x, y, z) -> T_init(z),
         S = S0)
    return nothing
end

# -------------------------
# State save helpers
# -------------------------
to_cpu(a) = Array(interior(a))

function save_state(filepath, model)
    jldopen(filepath, "w") do file
        file["time"] = model.clock.time
        file["iteration"] = model.clock.iteration
        file["u"] = to_cpu(model.velocities.u)
        file["v"] = to_cpu(model.velocities.v)
        file["w"] = to_cpu(model.velocities.w)
        file["T"] = to_cpu(model.tracers.T)
        file["S"] = to_cpu(model.tracers.S)
    end
    @info "Saved state to $filepath"
    return nothing
end

# -------------------------
# Progress
# -------------------------
function progress(sim)
    c = sim.model.clock
    @printf("iter: %d | time: %.4f yr | Δt: %.2f min\n",
            c.iteration, c.time / year, sim.Δt / minute)
    return nothing
end

# -------------------------
# Main run
# -------------------------
function run_simulation(; total_years::Float64 = 70.0,
                          Δt0 = 45minute,
                          archive_days::Int = 5,
                          prefix::String = "double_gyre_physics")

    @assert total_years > 1 "Need at least > 1 year so the final year can be archived."

    spinup_years = total_years - 1

    model = build_model()
    set_initial_conditions!(model)

    simulation = Simulation(model; Δt=Δt0, stop_time=spinup_years * year)

    wizard = TimeStepWizard(cfl=0.7,
                            max_change=1.1,
                            min_change=0.5,
                            max_Δt=60minute)

    simulation.callbacks[:wizard]   = Callback(wizard, IterationInterval(10))
    simulation.callbacks[:progress] = Callback(progress, IterationInterval(200))

    simulation.output_writers[:checkpointer] = Checkpointer(
        model;
        schedule = TimeInterval(30day),
        prefix = "$(prefix)_checkpoint",
        cleanup = true
    )

    @info "Phase 1: physics-only spin-up for $(spinup_years) years"
    run!(simulation)

    save_state("$(prefix)_cycle_init.jld2", model)

    simulation.stop_time = total_years * year

    simulation.output_writers[:offline_year] = JLD2Writer(
        model,
        (
            u = model.velocities.u,
            v = model.velocities.v,
            w = model.velocities.w,
            T = model.tracers.T,
            S = model.tracers.S,
        );
        filename = "$(prefix)_offline_year.jld2",
        schedule = TimeInterval(archive_days * day),
        overwrite_existing = true
    )

    @info "Phase 2: archiving final year every $(archive_days) days"
    run!(simulation; checkpoint_at_end=true)

    save_state("$(prefix)_terminal_state.jld2", model)

    return nothing
end

# -------------------------
# CLI
# -------------------------
total_years  = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 70.0
Δt0_minutes  = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 45.0
archive_days = length(ARGS) >= 3 ? parse(Int, ARGS[3])     : 5

run_simulation(; total_years = total_years,
               Δt0 = Δt0_minutes * minute,
               archive_days = archive_days)