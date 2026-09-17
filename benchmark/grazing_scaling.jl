using Agate
using BenchmarkTools
using OceanBioME: Biogeochemistry, PrescribedPhotosyntheticallyActiveRadiation
using Oceananigans
using Oceananigans.Architectures: CPU, GPU
using Oceananigans.Fields: FunctionField
using Oceananigans.Grids: Periodic, Bounded

const NiPiZD = Agate.Models.NiPiZD
const USE_GPU = "--gpu" in ARGS
const QUICK = "--quick" in ARGS
const PFT_COUNTS = QUICK ? (2, 4) : (2, 4, 8, 16)
const ARCH = USE_GPU ? GPU() : CPU()
const FLOAT = Float32

if USE_GPU
    @eval using CUDA
    CUDA.functional() || error("CUDA is not functional on this machine")
end

@inline benchmark_PAR(x, y, z, t) = FLOAT(100)

function matched_size_structure(n)
    phytoplankton = exp.(range(log(FLOAT(1)), log(FLOAT(20)); length=n))
    zooplankton = FLOAT(10) .* phytoplankton
    return (;
        phytoplankton=(P=phytoplankton,),
        zooplankton=(Z=zooplankton,),
    )
end

function benchmark_grid()
    # Large enough that GPU measurements are not purely launch latency, while still
    # keeping the benchmark light enough for routine developer/CI use.
    size = USE_GPU ? (32, 32, 16) : (8, 8, 8)
    return RectilinearGrid(
        ARCH, FLOAT;
        topology=(Periodic, Periodic, Bounded),
        size,
        x=(FLOAT(0), FLOAT(size[1])),
        y=(FLOAT(0), FLOAT(size[2])),
        z=(FLOAT(-size[3]), FLOAT(0)),
    )
end

function build_model(n)
    grid = benchmark_grid()
    clock = Clock(; time=zero(grid))
    PAR = FunctionField{Center,Center,Center}(benchmark_PAR, grid; clock)
    light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(PAR)
    bgc = NiPiZD.construct(; grid, size_structure=matched_size_structure(n))
    biogeochemistry = Biogeochemistry(bgc; light_attenuation)

    # Disable resolved advection so the scaling signal is dominated by tracer/BGC work.
    model = NonhydrostaticModel(
        grid;
        clock,
        biogeochemistry,
        advection=nothing,
    )

    plankton = merge(
        NamedTuple{Tuple(Symbol(:P_, i) for i in 1:n)}(ntuple(_ -> FLOAT(0.02), n)),
        NamedTuple{Tuple(Symbol(:Z_, i) for i in 1:n)}(ntuple(_ -> FLOAT(0.05), n)),
    )
    set!(model; N=FLOAT(7), D=FLOAT(0.01), plankton...)
    return model
end

sync_backend() = USE_GPU ? CUDA.synchronize() : nothing

function timed_step!(model)
    time_step!(model, FLOAT(60))
    sync_backend()
    return nothing
end

function benchmark_case(n)
    model = build_model(n)

    # Compile/warm once outside the measurement.
    timed_step!(model)

    trial = @benchmark timed_step!($model) samples=(QUICK ? 5 : 10) evals=1
    return BenchmarkTools.median(trial).time * 1e-9
end

backend_name = USE_GPU ? "GPU" : "CPU"
println("Agate NiPiZD grazing/PFT scaling benchmark ($backend_name)")
println("grid = ", size(benchmark_grid()), ", Float32, complete time_step!, advection=nothing")
println("P/Z\tedges\tmedian_ms\trelative_to_2P2Z")

base = nothing
for n in PFT_COUNTS
    seconds = benchmark_case(n)
    base === nothing && (base = seconds)
    println(n, "/", n, '\t', n * n, '\t', round(1e3 * seconds; digits=4), '\t', round(seconds / base; digits=3))
end
