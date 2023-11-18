using Agate
using Test

@testset "Agate.jl" begin
    @test Agate.placeholder_message() == "Hello, this is a placeholder for Agate.jl!"
    @test Agate.placeholder_message() != "Hello world!"
end
