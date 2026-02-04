"""Community helpers for common model families."""

"""Build a plankton community from a base community spec.

This helper updates only *structural* fields (`n` and `diameters`) for one or
more plankton groups, leaving all other group fields intact (e.g. `pft`).

The ordering of groups in the returned `NamedTuple` is the ordering of
`base` (i.e. `keys(base)`), which is stable and explicit.

Keywords
--------
- `n=(; ...)`: optional per-group size-class counts, keyed by group symbol.
- `diameters=(; ...)`: optional per-group diameter specifications, keyed by group symbol.

Examples
--------
```julia
community = build_plankton_community(base;
    n=(Z=10, P=20),
    diameters=(Z=zoo_diameters, P=phyto_diameters),
)
```
"""
function build_plankton_community(
    base::NamedTuple;
    n::NamedTuple=NamedTuple(),
    diameters::NamedTuple=NamedTuple(),
)
    base_keys = keys(base)

    # Explicit errors for unknown update keys.
    for k in keys(n)
        k in base_keys || throw(ArgumentError("n: unknown group symbol $k"))
    end
    for k in keys(diameters)
        k in base_keys || throw(ArgumentError("diameters: unknown group symbol $k"))
    end

    values_ = ntuple(length(base_keys)) do i
        g = base_keys[i]
        spec = getfield(base, g)

        # Only touch structural fields; preserve everything else.
        new_d = hasproperty(diameters, g) ? getfield(diameters, g) : getproperty(spec, :diameters)

        # `n` is optional when diameters are given explicitly as a vector/list.
        # Prefer an explicit `n` override, otherwise infer from explicit diameters,
        # otherwise fall back to the base spec.
        new_n = if hasproperty(n, g)
            getfield(n, g)
        elseif new_d isa AbstractVector
            length(new_d)
        elseif hasproperty(spec, :n)
            getproperty(spec, :n)
        else
            throw(
                ArgumentError(
                    "group $g: missing `n` and diameters are not an explicit vector; provide `n` explicitly",
                ),
            )
        end

        return (; spec..., n=new_n, diameters=new_d)
    end

    return NamedTuple{base_keys}(values_)
end
