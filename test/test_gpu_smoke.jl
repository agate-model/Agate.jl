using Agate
using Adapt
using Test

@testset "GPU smoke" begin
    cuda_ok = false
    try
        @eval using CUDA
        cuda_ok = CUDA.functional()
    catch
        cuda_ok = false
    end

    if !cuda_ok
        @test true
        return
    end

    bgc_type = Agate.Models.NiPiZD.construct(; FT=Float32)
    bgc = bgc_type()
    bgc_gpu = Adapt.adapt(CUDA.CuArray, bgc)

    out = CUDA.zeros(Float32, 1)

    function kernel!(out, bgc)
        i = CUDA.threadIdx().x
        if i == 1
            N = 7.0f0
            D = 1.0f0
            Z1 = 0.05f0
            Z2 = 0.05f0
            P1 = 0.01f0
            P2 = 0.01f0
            PAR = 100.0f0

            out[1] = bgc(Val(:N), 0.0f0, 0.0f0, 0.0f0, 0.0f0, N, D, Z1, Z2, P1, P2, PAR)
        end
        return
    end

    CUDA.@cuda threads=1 kernel!(out, bgc_gpu)
    @test isfinite(Array(out)[1])
end
