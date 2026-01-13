"""Expression-building helpers.

This submodule exists to avoid a circular dependency between `Agate.Utils` and
`Agate.Library`:

- `Agate.Utils` depends on many `Agate.Library.*` modules.
- Several `Agate.Library.*` modules need small `Expr` helpers (like `sum_expr`).

Keeping these helpers in `Agate.Library.ExprUtils` allows both `Library` modules
and `Utils` to import them without depending on each other.
"""

module ExprUtils

export sum_expr

"""
    sum_expr(terms)

Return an expression that sums a collection of AST terms.

Terms may be `Expr`, `Symbol`, or literal values (e.g. numbers).
If `terms` is empty, returns `0`.
"""
function sum_expr(terms::AbstractVector)
    isempty(terms) && return 0

    s = terms[1]
    for i in 2:length(terms)
        s = :($s + $(terms[i]))
    end
    return s
end

end # module ExprUtils
