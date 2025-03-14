using Test

@testset "Test Pauli Noise" begin
    nq = 8
    nl = 4
    W = Inf
    min_abs_coeff = 0.0

    opind = rand(1:nq)
    pstr = PauliString(nq, :Z, opind)

    topo = bricklayertopology(nq; periodic=false)
    circ = hardwareefficientcircuit(nq, nl; topology=topo)

    m = countparameters(circ)

    depolarizing_circ = deepcopy(circ)
    pauli_circ = deepcopy(circ)

    where_ind = rand(1:m)
    q_ind = opind
    noise_p = rand() * 0.2
    insert!(depolarizing_circ, where_ind, DepolarizingNoise(q_ind))
    insert!(pauli_circ, where_ind, PauliPropagation.PauliZDamping(q_ind))
    insert!(pauli_circ, where_ind, PauliPropagation.PauliYDamping(q_ind))
    insert!(pauli_circ, where_ind, PauliPropagation.PauliXDamping(q_ind))

    Random.seed!(42)
    thetas1 = randn(m)
    thetas2 = deepcopy(thetas1)
    insert!(thetas1, where_ind, noise_p)
    insert!(thetas2, where_ind, noise_p)
    insert!(thetas2, where_ind, noise_p)
    insert!(thetas2, where_ind, noise_p)

    dnum1 = propagate(depolarizing_circ, pstr, thetas1; max_weight=W, min_abs_coeff=min_abs_coeff)

    dnum2 = propagate(pauli_circ, pstr, thetas2; max_weight=W, min_abs_coeff=min_abs_coeff)

    @test overlapwithzero(dnum1) ≈ overlapwithzero(dnum2)
end


@testset "Test Dephasing Noise" begin
    nq = 8
    nl = 4
    W = Inf
    min_abs_coeff = 0.0

    opind = rand(1:nq)
    pstr = PauliString(nq, :Z, opind)

    topo = bricklayertopology(nq; periodic=false)
    circ = hardwareefficientcircuit(nq, nl; topology=topo)

    m = countparameters(circ)

    dephasing_circ = deepcopy(circ)
    pauli_circ = deepcopy(circ)

    where_ind = rand(1:m)
    q_ind = opind
    noise_p = rand() * 0.2
    insert!(dephasing_circ, where_ind, DephasingNoise(q_ind))
    insert!(pauli_circ, where_ind, PauliPropagation.PauliYDamping(q_ind))
    insert!(pauli_circ, where_ind, PauliPropagation.PauliXDamping(q_ind))

    Random.seed!(42)
    thetas1 = randn(m)
    thetas2 = deepcopy(thetas1)
    insert!(thetas1, where_ind, noise_p)
    insert!(thetas2, where_ind, noise_p)
    insert!(thetas2, where_ind, noise_p)

    dnum1 = propagate(dephasing_circ, pstr, thetas1; max_weight=W, min_abs_coeff=min_abs_coeff)

    dnum2 = propagate(pauli_circ, pstr, thetas2; max_weight=W, min_abs_coeff=min_abs_coeff)

    @test overlapwithzero(dnum1) ≈ overlapwithzero(dnum2)
end