using Agate

bgc = Agate.Models.NiPiZD.construct()
args = (0, 0, 0, 0, 7.0, 0.0, 0.05, 0.05, 0.01, 0.01, 100.0)
tendency() = bgc(Val(:P_1), args...)

for _ in 1:100
    tendency()
end

allocations = @allocated tendency()
iterations = 100_000
elapsed = @elapsed begin
    value = zero(tendency())
    for _ in 1:iterations
        value += tendency()
    end
    value
end

println("allocations_per_call=", allocations)
println("seconds_per_call=", elapsed / iterations)
