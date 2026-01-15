using Agate
using Test

using Agate.Library.Equations: ParamVar, sum_over

@testset "Equations: sum_over" begin
    # Basic symbolic sum of plain Symbols.
    items = [:x, :y, :z]
    ae = sum_over(items) do sym, _
        sym
    end
    @test ae.node == :((x + y) + z)

    # Requirements are accumulated when summing ParamVar references.
    p = ParamVar{:p}()
    ae2 = sum_over(1:2) do _, i
        p[i]
    end
    @test ae2.req.vectors == [:p]
    @test isempty(ae2.req.scalars)
    @test isempty(ae2.req.matrices)
end
