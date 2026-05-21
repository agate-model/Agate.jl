using SciMLBase: ODEProblem
using ..Introspection: tracer_names
using ..Library.Light: CyclicalPAR

export ode_problem
export ode_initial_state

"""
    ode_initial_state(initial, tracers) -> Vector

Create a state vector ordered like `tracers`.

`initial` may be an `AbstractVector`, a `NamedTuple`, an `AbstractDict`, or any
object with properties matching the tracer names.
"""
function ode_initial_state(initial::AbstractVector, tracers)
    length(initial) == length(tracers) || throw(DimensionMismatch(
        "initial vector has length $(length(initial)), expected $(length(tracers)) for tracers $tracers"
    ))
    return collect(initial)
end

function ode_initial_state(initial, tracers)
    return [_initial_value(initial, tracer) for tracer in tracers]
end

function _initial_value(initial::NamedTuple, tracer::Symbol)
    haskey(initial, tracer) || throw(ArgumentError("Missing initial condition for tracer $tracer."))
    return getproperty(initial, tracer)
end

function _initial_value(initial::AbstractDict, tracer::Symbol)
    if haskey(initial, tracer)
        return initial[tracer]
    elseif haskey(initial, string(tracer))
        return initial[string(tracer)]
    else
        throw(ArgumentError("Missing initial condition for tracer $tracer."))
    end
end

function _initial_value(initial, tracer::Symbol)
    hasproperty(initial, tracer) || throw(ArgumentError("Missing initial condition for tracer $tracer."))
    return getproperty(initial, tracer)
end

function _bgc_from_parameter(default_bgc, p)
    p === nothing && return default_bgc

    if p isa NamedTuple && haskey(p, :bgc)
        return p.bgc
    elseif p isa AbstractDict && haskey(p, :bgc)
        return p[:bgc]
    else
        return p
    end
end

"""
    ode_problem(bgc; initial, tspan, tracers=tracer_names(bgc), light, p=bgc,
                return_metadata=false)

Create a `SciMLBase.ODEProblem` for a constructed Agate biogeochemistry object.

The returned problem can be solved directly with OrdinaryDiffEq, remade inside a
SciML `EnsembleProblem`, or remade inside a Turing model.

# Keywords

- `initial`: Initial conditions as a `NamedTuple`, `Dict`, vector, or object with
  tracer-named properties.
- `tspan`: ODE time span, for example `(0.0, 365days)`.
- `tracers`: State-vector order. Defaults to `Agate.Introspection.tracer_names(bgc)`.
- `light`: Callable `(u, p, t) -> PAR`. Defaults to `CyclicalPAR(-10)(t)`.
- `p`: The parameter object attached to the `ODEProblem`. By default this is
  `bgc`, so `remake(problem; p=new_bgc)` swaps in a new constructed model for
  parameter sweeps or Turing.jl.
- `return_metadata`: When `true`, return `(problem=problem, tracers=tracers, u0=u0)`.

# Example

```julia
bgc = Agate.Models.NiPiZD.construct()
initial = (N=7.0, D=0.01, Z1=0.01, Z2=0.01, P1=0.01, P2=0.1)
problem = Agate.Diagnostics.ode_problem(bgc; initial, tspan=(0.0, 365days))
```
"""
function ode_problem(bgc; initial, tspan,
                     tracers=Tuple(tracer_names(bgc)),
                     light=(u, p, t) -> CyclicalPAR(-10)(t),
                     p=bgc,
                     return_metadata::Bool=false)
    tracers = Tuple(tracers)
    u0 = ode_initial_state(initial, tracers)

    function rhs!(du, u, p, t)
        active_bgc = _bgc_from_parameter(bgc, p)
        PAR = light(u, p, t)

        for (i, tracer) in enumerate(tracers)
            du[i] = active_bgc(Val(tracer), 0, 0, 0, t, u..., PAR)
        end

        return nothing
    end

    problem = ODEProblem(rhs!, u0, tspan, p)
    return return_metadata ? (; problem, tracers, u0) : problem
end
