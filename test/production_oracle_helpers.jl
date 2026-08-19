using Agate
using Agate.Library.Light: FunctionFieldPAR
using OceanBioME: Biogeochemistry, BoxModel, BoxModelGrid
using Oceananigans: set!, time_step!
using Oceananigans.Units: day, minutes

const NIPIZD_PRODUCTION_TRACERS = (:N, :D, :Z_1, :Z_2, :P_1, :P_2)

function build_nipizd_production_box(
    ::Type{T}=Float64; sinking_tracers=nothing, open_bottom::Bool=true
) where {T<:AbstractFloat}
    grid = BoxModelGrid(T)
    bgc = Agate.Models.NiPiZD.construct(; grid, sinking_tracers, open_bottom)
    biogeochemistry = Biogeochemistry(bgc; light_attenuation=FunctionFieldPAR(; grid))
    return BoxModel(; biogeochemistry, grid), bgc
end

function set_nipizd_production_state!(box_model, ::Type{T}=Float64; ulp_perturb=false) where {T<:AbstractFloat}
    p1 = T(0.01)
    ulp_perturb && (p1 = nextfloat(p1))
    set!(
        box_model;
        N=T(7),
        D=zero(T),
        Z_1=T(0.05),
        Z_2=T(0.05),
        P_1=p1,
        P_2=T(0.01),
    )
    return box_model
end

function sample_box_fields(box_model, tracers::Tuple)
    return NamedTuple{tracers}(
        ntuple(length(tracers)) do i
            tracer = tracers[i]
            box_model.fields[tracer].data[1, 1, 1]
        end,
    )
end

function production_trajectory(
    box_model,
    tracers::Tuple;
    dt=10minutes,
    nsteps=Int(30day ÷ dt),
    sample_every=Int(12 * 60minutes ÷ dt),
)
    rows = NamedTuple[]
    push!(rows, (; step=0, time=0.0, sample_box_fields(box_model, tracers)...))

    for step in 1:nsteps
        time_step!(box_model, dt)
        if step % sample_every == 0 || step == nsteps
            push!(
                rows,
                (; step, time=Float64(step * dt), sample_box_fields(box_model, tracers)...),
            )
        end
    end
    return rows
end

function max_trajectory_difference(reference, perturbed, tracers::Tuple)
    length(reference) == length(perturbed) || throw(ArgumentError("trajectory lengths differ"))
    max_abs = 0.0
    max_rel = 0.0
    for tracer in tracers
        scale = maximum(
            max(abs(Float64(getproperty(reference[i], tracer))), abs(Float64(getproperty(perturbed[i], tracer))))
            for i in eachindex(reference, perturbed)
        )
        scale = max(scale, eps(Float64))
        for i in eachindex(reference, perturbed)
            a = Float64(getproperty(reference[i], tracer))
            b = Float64(getproperty(perturbed[i], tracer))
            delta = abs(a - b)
            max_abs = max(max_abs, delta)
            max_rel = max(max_rel, delta / scale)
        end
    end
    return (; max_abs, max_rel)
end

function write_trajectory_reference(path, rows, tracers::Tuple; metadata=(;))
    mkpath(dirname(path))
    open(path, "w") do io
        for (key, value) in pairs(metadata)
            println(io, "# ", key, "=", value)
        end
        println(io, join(("step", "time", string.(tracers)...), ','))
        for row in rows
            values = (row.step, row.time, (getproperty(row, tracer) for tracer in tracers)...)
            println(io, join(values, ','))
        end
    end
    return path
end

function read_trajectory_reference(path, tracers::Tuple)
    metadata = Dict{Symbol,Float64}()
    rows = NamedTuple[]
    header_seen = false
    for line in eachline(path)
        isempty(line) && continue
        if startswith(line, "# ")
            key, value = split(line[3:end], '='; limit=2)
            metadata[Symbol(key)] = parse(Float64, value)
            continue
        end
        if !header_seen
            expected = ["step", "time", string.(tracers)...]
            split(line, ',') == expected || throw(ArgumentError("unexpected trajectory header in $path"))
            header_seen = true
            continue
        end
        values = split(line, ',')
        step = parse(Int, values[1])
        time = parse(Float64, values[2])
        tracer_values = ntuple(length(tracers)) do i
            parse(Float64, values[i + 2])
        end
        push!(rows, (; step, time, NamedTuple{tracers}(tracer_values)...))
    end
    header_seen || throw(ArgumentError("trajectory reference has no header: $path"))
    return (; metadata, rows)
end
