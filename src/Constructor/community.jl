"""Community helpers for common model families."""

"""Build a `(Z, P)` community from a base community spec.

This helper is shared by NiPiZD and DARWIN, which both use the `(Z, P)` group set.

It only updates structural fields (`n` and `diameters`) and leaves the rest of the
group specification intact (e.g. `pft`).
"""
function build_ZP_community(
    base::NamedTuple;
    n_zoo::Int,
    n_phyto::Int,
    zoo_diameters,
    phyto_diameters,
)
    Z = (; base.Z..., n=n_zoo, diameters=zoo_diameters)
    P = (; base.P..., n=n_phyto, diameters=phyto_diameters)
    return (Z=Z, P=P)
end
