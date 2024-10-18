function paulinoise(nq, nl, W, min_abs_coeff)

    op = PauliString(nq, :Z, round(Int, nq / 2))

    topo = bricklayertopology(nq; periodic=false)
    circ = hardwareefficientcircuit(nq, nl; topology=topo)

    fastcirc = tofastgates(circ)
    m = countparameters(fastcirc)

    depolarizing_circ = deepcopy(fastcirc)
    pauli_circ = deepcopy(fastcirc)

    where_ind = rand(1:m)
    q_ind = rand(1:nq)
    noise_p = randn() * 0.2
    insert!(depolarizing_circ, where_ind, DepolarizingNoise(q_ind))
    insert!(pauli_circ, where_ind, PauliZNoise(q_ind))
    insert!(pauli_circ, where_ind, PauliYNoise(q_ind))
    insert!(pauli_circ, where_ind, PauliXNoise(q_ind))

    Random.seed!(42)
    thetas1 = randn(m)
    thetas2 = deepcopy(thetas1)
    insert!(thetas1, where_ind, noise_p)
    insert!(thetas2, where_ind, noise_p)
    insert!(thetas2, where_ind, noise_p)
    insert!(thetas2, where_ind, noise_p)

    dnum1 = mergingbfs(depolarizing_circ, op, thetas1; max_weight=W, min_abs_coeff=min_abs_coeff)

    dnum2 = mergingbfs(pauli_circ, op, thetas2; max_weight=W, min_abs_coeff=min_abs_coeff)

    return overlapwithzero(dnum1) ≈ overlapwithzero(dnum2)
end