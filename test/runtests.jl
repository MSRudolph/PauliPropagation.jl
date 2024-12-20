using PauliPropagation
using Test
using Random

@testset "PauliPropagation.jl" begin

    include("test_propagate.jl")
    @test numericalPP(8, 4, Inf, 0.0) ≈ 0.21720058439757214
    @test hybridPP(8, 4, Inf, 0.0, Inf) ≈ 0.21720058439757214
    @test surrogatePP(8, 4, Inf, Inf) ≈ 0.21720058439757214

    include("test_datatypes.jl")
    @test isa(createpaulistring(7), PauliString)
    @test isa(createpaulisum(21), PauliSum)
    @test isa(addtopaulisum(65), PauliSum)

    include("test_paulialgebra_utils.jl")

    Random.seed!(42)
    include("test_noisechannels.jl")
    @test paulinoise(8, 4, Inf, 0.0)
    @test dephasingnoise(17, 3, Inf, 0.0)

    include("test_cliffordgates.jl")

    include("test_frozengates.jl")

    include("test_overlaps.jl")

    include("test_paulirotations.jl")

    include("test_paulioperations.jl")

    include("test_truncations.jl")

    include("test_numericalcertificates.jl")

end
