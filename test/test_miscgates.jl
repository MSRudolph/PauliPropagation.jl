using LinearAlgebra
using Random
using Test

@testset "Unitaries for Miscillaneous Gates Tests" begin
    """Test the unitary matrices."""
    tol = 1e-12

    # Test using single-qubit PauliRotation gate
    @testset "PauliRotation Y" begin
        pauligate = PauliRotation(:Y, 1)
        theta = Random.randn()

        U = tomatrix(pauligate, theta)
        exptected_U = [
            [cos(theta / 2) -sin(theta / 2)];
            [sin(theta / 2) cos(theta / 2)]
        ]
        @test LinearAlgebra.norm(U - exptected_U) < tol
    end

    # Test using T gate
    @testset "TUnitary" begin
        tgate = PauliPropagation.TGate(1)
        U = tomatrix(tgate)
        expected_U = [
            [1 0];
            [0 exp(1.0im * pi / 4)]
        ]
        @test LinearAlgebra.norm(U - expected_U) < tol
    end
end
