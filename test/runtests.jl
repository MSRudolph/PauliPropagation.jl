using PauliPropagation
using Test

@testset "PauliPropagation.jl" begin
    # Write your tests here.
    include("test_apply.jl")
    th = randn()
    @test all(applypauligate([:Z], [1], th) .== (0x0003, 1.0))
    th = randn()
    @test all(applypauligate([:X], [1], th) .≈ (0x0003, cos(th), 0x0002, sin(th)))
    th = randn()
    @test all(applypauligate([:Y], [1], th) .≈ (0x0003, cos(th), 0x0001, -sin(th)))


    include("test_mergingbfs.jl")
    @test numericalPP(8, 4, Inf, 0.0) ≈ 0.21720058439757214
    @test hybridPP(8, 4, Inf, 0.0, Inf) ≈ 0.21720058439757214
    @test surrogatePP(8, 4, Inf, Inf) ≈ 0.21720058439757214


    Random.seed!(42)
    include("test_noisechannels.jl")
    @test paulinoise(8, 4, Inf, 0.0)

end
