using Test

@testset "Test FrozenGates" begin

    nq = 9

    nl = 3

    circ = efficientsu2circuit(nq, nl)

    insert!(circ, 6, PauliPropagation.PauliXDamping(1))
    insert!(circ, 10, DepolarizingNoise(3))
    insert!(circ, 21, AmplitudeDampingNoise(7))

    m = countparameters(circ)

    PauliXDamping_index = getparameterindices(circ, PauliPropagation.PauliXDamping)
    depolarizingnoise_index = getparameterindices(circ, DepolarizingNoise)
    amplitudedampingnoise_index = getparameterindices(circ, AmplitudeDampingNoise)

    thetas = randn(m)
    thetas[PauliXDamping_index] .= 0.12
    thetas[depolarizingnoise_index] .= 0.23
    thetas[amplitudedampingnoise_index] .= 0.02


    frozen_circ = freeze(circ, thetas)

    pstr = PauliString(nq, [:Y, :X], [3, 6])

    psum1 = propagate(circ, pstr, thetas)
    psum2 = propagate(frozen_circ, pstr)

    @test overlapwithzero(psum1) ≈ overlapwithzero(psum2)

end
