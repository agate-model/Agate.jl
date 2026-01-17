module Mortality

using ...Equations: AExpr, req, _to_aexpr, merge_requirements, sum_over

export linear_loss, quadratic_loss, linear_loss_sum, quadratic_loss_sum

@inline _vec(PV, key::Symbol, idx::Int) = getproperty(PV, key)[idx]

"""
    linear_loss(PV, plankton_sym, idx; mortality=:linear_mortality)

Construction-time linear loss term `PV[mortality][idx] * plankton_sym`.
"""
@inline function linear_loss(
    PV,
    plankton_sym::Symbol,
    idx::Int;
    mortality::Symbol=:linear_mortality,
)
    return _vec(PV, mortality, idx) * plankton_sym
end

"""
    quadratic_loss(PV, plankton_sym, idx; mortality=:quadratic_mortality)

Construction-time quadratic loss term `PV[mortality][idx] * plankton_sym^2`.
"""
@inline function quadratic_loss(
    PV,
    plankton_sym::Symbol,
    idx::Int;
    mortality::Symbol=:quadratic_mortality,
)
    return _vec(PV, mortality, idx) * plankton_sym * plankton_sym
end

"""
    linear_loss_sum(PV, plankton_syms; mortality=:linear_mortality)

Sum linear losses over a set of plankton tracer symbols.
"""
function linear_loss_sum(
    PV,
    plankton_syms::AbstractVector{Symbol};
    mortality::Symbol=:linear_mortality,
)
    return sum_over(plankton_syms) do sym, i
        _vec(PV, mortality, i) * sym
    end
end

"""
    quadratic_loss_sum(PV, plankton_syms; mortality=:quadratic_mortality)

Sum quadratic losses over a set of plankton tracer symbols.
"""
function quadratic_loss_sum(
    PV,
    plankton_syms::AbstractVector{Symbol};
    mortality::Symbol=:quadratic_mortality,
)
    return sum_over(plankton_syms) do sym, i
        _vec(PV, mortality, i) * sym * sym
    end
end

end # module
