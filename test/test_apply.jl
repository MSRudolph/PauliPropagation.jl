# TODO: Should be test gates
# YT: these have been moved to test/test_pauligates.jl
function applypauligate(pauli_generator, qinds, theta)
    nq = 5
    symbs = [:I for _ in 1:nq]
    symbs[1] = :Z   # as symbol. Also works but is slower.
    pstr = symboltoint(symbs)

    gate = PauliGate(pauli_generator, qinds)
    fastgate = tofastgates(gate, nq)
    @test apply(gate, pstr, theta) == apply(fastgate, pstr, theta)
    return apply(fastgate, pstr, theta)
end