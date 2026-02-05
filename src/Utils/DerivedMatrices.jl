"""Utilities for *derived* interaction matrices.

Some models expose low-level interaction *traits* (e.g. predator-prey size ratio
preferences) and derive interaction matrices (palatability, assimilation) from
those traits.

Derived matrices are recomputed during construction when:

- a matrix is *not* explicitly overridden, and
- at least one of its declared dependencies *is* explicitly overridden.

This allows users to tweak traits (small surface, easy to reason about) without
needing to hand-build full matrices.
"""

"""A derived-matrix provider and the parameter keys it depends on."""
struct MatrixProvider{F,Deps}
    f::F
    deps::Deps
    function MatrixProvider(f; deps=())
        deps_norm = if deps === nothing
            ()
        elseif deps isa Symbol
            (deps,)
        elseif deps isa AbstractVector
            Tuple(deps)
        else
            deps
        end
        return new{typeof(f),typeof(deps_norm)}(f, deps_norm)
    end
end
"""Return a `NamedTuple` mapping matrix keys to `MatrixProvider`s.

Factories override this to declare which matrices can be derived from other
parameters.
"""
derived_matrix_specs(::AbstractBGCFactory) = (;)

"""Resolve derived matrices into `params`.

`explicit_override_keys` should include keys explicitly provided by the user
either through `parameters` or `interactions`.

If a derived matrix is *explicitly* overridden, its value is kept as-is.
Otherwise, the matrix is recomputed if any of its dependencies are explicitly
overridden.
"""
function resolve_derived_matrices(
    factory::AbstractBGCFactory,
    context::CommunityContext,
    params::NamedTuple,
    explicit_override_keys::Tuple{Vararg{Symbol}},
)
    derived_matrix_providers = derived_matrix_specs(factory)
    isempty(derived_matrix_providers) && return params

    override_set = Set(explicit_override_keys)
    matrix_keys = propertynames(derived_matrix_providers)

    vals = ntuple(
        i -> begin
            key = matrix_keys[i]
            provider = getproperty(derived_matrix_providers, key)

            # If the matrix itself was explicitly overridden, keep it.
            if key in override_set && hasproperty(params, key)
                return getproperty(params, key)
            end

            # Only recompute when a declared dependency was explicitly overridden.
            if !isempty(provider.deps) && !any(d -> d in override_set, provider.deps)
                return if hasproperty(params, key)
                    getproperty(params, key)
                else
                    provider.f(factory, context, params)
                end
            end

            return provider.f(factory, context, params)
        end,
        length(matrix_keys),
    )

    updates = NamedTuple{matrix_keys}(vals)
    return merge(params, updates)
end
