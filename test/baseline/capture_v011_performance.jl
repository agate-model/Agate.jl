using Agate
using JSON
using OceanBioME: BoxModelGrid
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers

const NiPiZD = Agate.Models.NiPiZD

GC.gc()
construction = @timed NiPiZD.construct(; grid=BoxModelGrid())
bgc = construction.value

tracers = required_biogeochemical_tracers(bgc)
state = (7.0, 1.0, 0.05, 0.05, 0.01, 0.01)
PAR = 100.0

GC.gc()
first_tendency = @timed bgc(Val(:P_1), 0, 0, 0, 0, state..., PAR)

result = Dict(
    "package" => Dict(
        "name" => "Agate",
        "version" => string(Base.pkgversion(Agate)),
    ),
    "model" => "NiPiZD",
    "julia_version" => string(VERSION),
    "cpu" => Sys.CPU_NAME,
    "threads" => Threads.nthreads(),
    "scalar_type" => "Float64",
    "tracer_order" => String.(tracers),
    "construction" => Dict(
        "seconds" => construction.time,
        "bytes" => construction.bytes,
    ),
    "first_tendency_call" => Dict(
        "tracer" => "P_1",
        "seconds" => first_tendency.time,
        "bytes" => first_tendency.bytes,
        "value" => first_tendency.value,
    ),
)

print(JSON.json(result, 2))
println()
