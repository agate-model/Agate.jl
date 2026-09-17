using Agate
using OceanBioME:
    Biogeochemistry, BoxModelGrid, BoxModel, PrescribedPhotosyntheticallyActiveRadiation
using Oceananigans: set!, time_step!
using Oceananigans.Fields: ConstantField
using Oceananigans.Units: day
using Test

const NiPiZDReference = Agate.Models.NiPiZD
const NIPIZD_REFERENCE_TRACERS = (:N, :D, :Z_1, :Z_2, :P_1, :P_2)
const NIPIZD_REFERENCE_COLUMNS = (:time_days, :total_nitrogen, NIPIZD_REFERENCE_TRACERS...)
const NIPIZD_REFERENCE_INITIAL = (
    N=7.0,
    D=0.01,
    Z_1=0.01,
    Z_2=0.02,
    P_1=0.01,
    P_2=0.1,
)
const NIPIZD_REFERENCE_DT = day / 48
const NIPIZD_REFERENCE_SAVE_EVERY = 12
const NIPIZD_REFERENCE_DURATION = 60day
const NIPIZD_REFERENCE_RTOL = 1e-12
const NIPIZD_REFERENCE_ATOL = 1e-14

const NIPIZD_REFERENCE_PATH = joinpath(
    @__DIR__, "reference", "nipizd_v0.14.0_production_reference.csv"
)

function nipizd_reference_row(time_days, box_model)
    tracer_values = ntuple(length(NIPIZD_REFERENCE_TRACERS)) do i
        tracer = NIPIZD_REFERENCE_TRACERS[i]
        Float64(getproperty(box_model.fields, tracer).data[1, 1, 1])
    end
    state = NamedTuple{NIPIZD_REFERENCE_TRACERS}(tracer_values)
    return (;
        time_days=Float64(time_days),
        total_nitrogen=sum(tracer_values),
        state...,
    )
end

function simulate_nipizd_reference()
    grid = BoxModelGrid()
    bgc = NiPiZDReference.construct()
    light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(ConstantField(100.0))
    bgc_model = Biogeochemistry(bgc; light_attenuation)
    box_model = BoxModel(; biogeochemistry=bgc_model)
    set!(box_model; NIPIZD_REFERENCE_INITIAL...)

    sample_dt = NIPIZD_REFERENCE_SAVE_EVERY * NIPIZD_REFERENCE_DT
    n_samples = Int(NIPIZD_REFERENCE_DURATION / sample_dt)
    rows = NamedTuple[nipizd_reference_row(0.0, box_model)]

    for sample in 1:n_samples
        for _ in 1:NIPIZD_REFERENCE_SAVE_EVERY
            time_step!(box_model, NIPIZD_REFERENCE_DT)
        end
        push!(rows, nipizd_reference_row(sample * sample_dt / day, box_model))
    end

    return rows
end

function write_nipizd_reference(path=NIPIZD_REFERENCE_PATH)
    rows = simulate_nipizd_reference()
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(string.(NIPIZD_REFERENCE_COLUMNS), ","))
        for row in rows
            values = (repr(getproperty(row, name)) for name in NIPIZD_REFERENCE_COLUMNS)
            println(io, join(values, ","))
        end
    end
    return path
end

function read_nipizd_reference(path=NIPIZD_REFERENCE_PATH)
    isfile(path) || error(
        "NiPiZD reference trajectory is missing at $path. " *
        "Generate it with `julia --project=. test/test_nipizd_reference_trajectory.jl --generate`."
    )
    reference_rows = filter(
        row -> !isempty(row) && !startswith(row, '#'),
        readlines(path),
    )
    header = Symbol.(split(first(reference_rows), ','))
    reference = map(reference_rows[2:end]) do row
        values = parse.(Float64, split(row, ','))
        NamedTuple{Tuple(header)}(Tuple(values))
    end
    return header, reference
end

if "--generate" in ARGS
    path = write_nipizd_reference()
    @info "Wrote NiPiZD reference trajectory" path
else
    header, reference = read_nipizd_reference()

    @testset "NiPiZD v0.14 production trajectory" begin
        @test Tuple(header) == NIPIZD_REFERENCE_COLUMNS
        @test [row.time_days for row in reference] == collect(0.0:0.25:60.0)

        actual = simulate_nipizd_reference()
        @test length(actual) == length(reference)

        mismatches = NamedTuple[]
        for (actual_row, expected_row) in zip(actual, reference)
            for tracer in NIPIZD_REFERENCE_TRACERS
                actual_value = getproperty(actual_row, tracer)
                expected_value = getproperty(expected_row, tracer)
                isapprox(
                    actual_value,
                    expected_value;
                    rtol=NIPIZD_REFERENCE_RTOL,
                    atol=NIPIZD_REFERENCE_ATOL,
                ) && continue

                push!(
                    mismatches,
                    (;
                        time_days=expected_row.time_days,
                        tracer,
                        actual=actual_value,
                        expected=expected_value,
                        absolute_error=abs(actual_value - expected_value),
                    ),
                )
            end
        end

        if !isempty(mismatches)
            worst = mismatches[argmax(getproperty.(mismatches, :absolute_error))]
            @info "NiPiZD trajectory mismatch" count=length(mismatches) time_days=worst.time_days tracer=worst.tracer actual=worst.actual expected=worst.expected absolute_error=worst.absolute_error
        end

        @test isempty(mismatches)
    end
end
