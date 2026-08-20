using Agate
using Agate.Library.Light: FunctionFieldPAR
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers
using Oceananigans.Units: day, minutes

const NiPiZD = Agate.Models.NiPiZD

function generate_nipizd_production_reference()
    grid = BoxModelGrid()
    bgc = NiPiZD.construct(; grid)
    model = BoxModel(;
        biogeochemistry=Biogeochemistry(
            bgc; light_attenuation=FunctionFieldPAR(; grid)
        ),
        grid,
    )
    set!(model; N=7.0, P_1=0.01, P_2=0.01, Z_1=0.05, Z_2=0.05, D=0.0)

    tracers = required_biogeochemical_tracers(bgc)
    dt = 30minutes
    sample_interval = 1day
    stop_time = 180day
    sample_steps = Int(sample_interval / dt)
    nsteps = Int(stop_time / dt)

    version = Base.pkgversion(Agate)
    reference_dir = joinpath(@__DIR__, "reference")
    mkpath(reference_dir)
    path = joinpath(reference_dir, "nipizd_production_v$(version)_candidate.csv")

    open(path, "w") do io
        println(io, "# model = Agate.Models.NiPiZD.construct")
        println(io, "# version = $version")
        println(io, "# integration = BoxModel, dt=30 minutes, stop=180 days")
        println(io, "# sampling = every 1 day")
        println(io, "# initial_state = N=7.0, P_1=0.01, P_2=0.01, Z_1=0.05, Z_2=0.05, D=0.0")
        println(io, "# light = FunctionFieldPAR on BoxModelGrid")
        println(io, join(("day", string.(tracers)...), ','))

        for step in 1:nsteps
            time_step!(model, dt)
            step % sample_steps == 0 || continue
            values = Tuple(model.fields[name].data[1, 1, 1] for name in tracers)
            day_value = step * dt / day
            println(io, join((string(day_value), string.(values)...), ','))
        end
    end

    println(path)
    return path
end

generate_nipizd_production_reference()
