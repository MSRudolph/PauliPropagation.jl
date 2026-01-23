## Test heisenberg=true and heisenberg=false equivalence for implemented gates
using Random

@testset "Test propagate heisenberg vs schrodinger equivalence" begin
    nq = 5
    nl1 = 2
    nl2 = 3
    nl_half = max(nl1, nl2)

    min_abs_coeff = 0

    pstr_beginning = PauliString(nq, :Z, 2)
    pstr_end = PauliString(nq, [:X, :Y], [1, 3])

    topology = staircasetopology(nq; periodic=true)

    # contains Pauli rotations and CNOTs
    circ1 = hardwareefficientcircuit(nq, nl1; topology=topology)
    circ2 = tiltedtfitrottercircuit(nq, nl2; topology=topology)

    # tests frozen gates and noise
    noise_layer = Gate[DepolarizingNoise(i, 0.1) for i in 1:nq]

    circ = vcat(circ1, noise_layer, circ2)

    thetas1 = randn(countparameters(circ1))
    thetas2 = randn(countparameters(circ2))
    thetas = vcat(thetas1, thetas2)

    psum_beginning = propagate(circ, pstr_beginning, thetas; heisenberg=false, min_abs_coeff)
    val1 = scalarproduct(psum_beginning, pstr_end)

    psum_end = propagate!(circ, pstr_end, thetas; heisenberg=true, min_abs_coeff)
    val2 = scalarproduct(pstr_beginning, psum_end)

    psum_beginning = propagate(circ1, pstr_beginning, thetas1; heisenberg=false, min_abs_coeff)
    psum_beginning = propagate!(noise_layer, psum_beginning; heisenberg=false, min_abs_coeff)
    psum_end = propagate(circ2, pstr_end, thetas2; heisenberg=true, min_abs_coeff)
    val3 = scalarproduct(psum_beginning, psum_end)

    @test val1 ≈ val2 ≈ val3

end

@testset "Test inversion" begin
    nx = 3
    ny = 2
    nq = nx * ny
    nl = 3
    topology = rectangletopology(nx, ny; periodic=false)

    min_abs_coeff = 1e-12

    circ = hardwareefficientcircuit(nq, nl; topology=topology)
    thetas = randn(countparameters(circ))

    fcirc = freeze(circ, thetas)

    pstr = PauliString(nq, rand([:X, :Y, :Z], nq), 1:nq)
    psum_forward = propagate(fcirc, pstr; heisenberg=false, min_abs_coeff)
    psum_backward = propagate(fcirc, psum_forward; heisenberg=true, min_abs_coeff)
    @test PauliSum(pstr) ≈ psum_backward

    # flip direction
    pstr = PauliString(nq, rand([:X, :Y, :Z], nq), 1:nq)
    psum_forward = propagate(fcirc, pstr; heisenberg=true, min_abs_coeff)
    psum_backward = propagate(fcirc, psum_forward; heisenberg=false, min_abs_coeff)
    @test PauliSum(pstr) ≈ psum_backward

end
