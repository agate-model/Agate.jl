using Agate
using OceanBioME:
    Biogeochemistry, BoxModel, PrescribedPhotosyntheticallyActiveRadiation
using Oceananigans: set!, time_step!
using Oceananigans.Fields: ConstantField
using Oceananigans.Units: day

const tracers = (:N, :D, :Z_1, :Z_2, :P_1, :P_2)
const columns = (:time_days, :total_nitrogen, tracers...)
const initial = (N=7.0, D=0.01, Z_1=0.01, Z_2=0.02, P_1=0.01, P_2=0.1)
const dt = day / 48
const save_every = 12
const output = joinpath(
    @__DIR__, "..", "test", "reference", "nipizd_definition_v0.2.0_reference.csv"
)

function state_row(time_days, model)
    values = ntuple(length(tracers)) do i
        Float64(getproperty(model.fields, tracers[i]).data[1, 1, 1])
    end
    state = NamedTuple{tracers}(values)
    return (; time_days, total_nitrogen=sum(values), state...)
end

bgc, recipe = Agate.Models.NiPiZD.construct_plus_recipe()
light = PrescribedPhotosyntheticallyActiveRadiation(ConstantField(100.0))
model = BoxModel(; biogeochemistry=Biogeochemistry(bgc; light_attenuation=light))
set!(model; initial...)

rows = [state_row(0.0, model)]
for sample in 1:240
    for _ in 1:save_every
        time_step!(model, dt)
    end
    push!(rows, state_row(sample / 4, model))
end

mkpath(dirname(output))
open(output, "w") do io
    println(io, "# generator=scripts/generate_nipizd_reference.jl")
    println(io, "# agate_version=", Base.pkgversion(Agate))
    println(io, "# nipizd_definition_version=", recipe.definition_version)
    println(io, "# julia_version=", VERSION)
    println(io, "# constructor=Agate.Models.NiPiZD.construct_plus_recipe()")
    println(io, "# forcing=constant_PAR_100.0")
    println(io, "# dt_seconds=", dt)
    println(io, "# save_interval_seconds=", dt * save_every)
    initial_conditions = join(
        (string(name, "=", getproperty(initial, name)) for name in keys(initial)), ","
    )
    println(io, "# initial_conditions=", initial_conditions)
    println(io, "# tracer_order=", join(string.(tracers), ","))
    println(io, join(string.(columns), ","))
    for row in rows
        println(io, join((repr(getproperty(row, name)) for name in columns), ","))
    end
end

@info "Wrote NiPiZD reference trajectory" output
