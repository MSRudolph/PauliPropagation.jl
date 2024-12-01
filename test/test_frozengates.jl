using Test

@testset "Test FrozenGates" begin

    nq = 9

    nl = 3

    circ = efficientsu2circuit(nq, nl)

    insert!(circ, 6, PauliXNoise(1))
    insert!(circ, 10, DepolarizingNoise(3))
    insert!(circ, 21, AmplitudeDampingNoise(7))

    m = countparameters(circ)

    thetas = randn(m) .* 0.2

    frozen_circ = freeze(circ, thetas)

    op = PauliString(nq, [:Y, :X], [3, 6])

    psum1 = propagate(circ, op, thetas)
    psum2 = propagate(frozen_circ, op, [0])

    @test overlapwithzero(psum1) ≈ overlapwithzero(psum2)

end
