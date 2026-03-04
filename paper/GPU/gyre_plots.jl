using Oceananigans
using Oceananigans.Units: day
using Statistics: mean

const year = 365day
T = FieldTimeSeries("double_gyre_online_noCC_noFW.jld2", "T")
xC,yC,zC = nodes(T.grid, Center(), Center(), Center())
k = length(zC)

for n in max(1, length(T.times)-6):length(T.times)
    A = interior(T[n])[:, :, k]
    println("n=$n  t=$(round(T.times[n]/day,digits=1)) d  mean=$(mean(A))  min=$(minimum(A))  max=$(maximum(A))")
end