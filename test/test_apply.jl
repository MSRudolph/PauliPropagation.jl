# TODO: Should be test gates
# YT: these have been moved to test/test_pauligates.jl
function applypauligate(pauli_generator, qinds, theta)
    nq = 5
    symbs = [:I for _ in 1:nq]
    symbs[1] = :Z   # as symbol. Also works but is slower.
    op = symboltoint(symbs)

    gate = PauliGate(pauli_generator, qinds)
    fastgate = tofastgates(gate, nq)
    @test apply(gate, op, theta) == apply(fastgate, op, theta)
    return apply(fastgate, op, theta)
end