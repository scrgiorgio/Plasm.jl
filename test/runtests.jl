using Plasm
using Test

@testset "Plasm.jl" begin
    @test greet_your_package_name() == "Hello"
    @test greet_your_package_name() != "Hi"
end
