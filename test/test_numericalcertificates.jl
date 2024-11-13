@testset "numerical certificate Tests" begin
    nq = rand(1:100)

    nl = 4
    topo = bricklayertopology(nq; periodic=false)

    circ = efficientsu2circuit(nq, nl; topology=topo)
    nparams = countparameters(circ)
    ngates = length(circ)


    # Numerical Coefficient
    op = PauliString(nq, :X, rand(1:nq))
    @test isa(estimateaverageerror(circ, op, 100), Float64)
    op = PauliString(nq, :Y, rand(1:nq))
    @test isa(estimateaverageerror(circ, op, 100, 0.87236582), Float64)

    # Weight Truncation
    op = PauliString(nq, :Z, rand(1:nq))
    @test estimateaverageerror(circ, op, 100, π; max_weight=0) >= estimateaverageerror(circ, op, 100, π; max_weight=nq) == 0.0
    @test estimateaverageerror(circ, op, 100, ones(nparams) * 1.23; max_weight=0) >= estimateaverageerror(circ, op, 100, ones(nparams) * 1.23; max_weight=nq) == 0.0

    # Frequency Truncation
    wop = wrapcoefficients(PauliString(nq, :Z, rand(1:nq)))
    @test estimateaverageerror(circ, wop, 100, π; max_freq=0) >= estimateaverageerror(circ, op, 100, π; max_freq=nq) == 0.0
    @test estimateaverageerror(circ, wop, 100, ones(nparams) * 1.23; max_freq=0) >= estimateaverageerror(circ, op, 100, ones(nparams) * 1.23; max_freq=nq) == 0.0

    # Small-angle Truncation
    wop = wrapcoefficients(PauliString(nq, :Z, rand(1:nq)))
    @test estimateaverageerror(circ, wop, 100, π; max_sins=0) >= estimateaverageerror(circ, op, 100, π; max_sins=nq) == 0.0
    @test estimateaverageerror(circ, wop, 100, ones(nparams) * 1.23; max_sins=0) >= estimateaverageerror(circ, op, 100, ones(nparams) * 1.23; max_sins=nq) == 0.0

    # Numerical Coefficient
    op = PauliString(nq, :X, rand(1:nq))
    @test typeof(montecarlopropagation(circ, op)) <: Tuple{typeof(op),Bool}
    op = PauliString(nq, :Y, rand(1:nq))
    @test typeof(montecarlopropagation(circ, op, 0.5)) <: Tuple{typeof(op),Bool}
    op = PauliString(nq, :Z, rand(1:nq))
    @test typeof(montecarlopropagation(circ, op, ones(nparams) * π)) <: Tuple{typeof(op),Bool}

    # NumericPathProperties Coefficient
    wop = wrapcoefficients(PauliString(nq, :X, rand(1:nq)))
    @test typeof(montecarlopropagation(circ, wop)) <: Tuple{typeof(wop),Bool}
    wop = wrapcoefficients(PauliString(nq, :Y, rand(1:nq)))
    @test typeof(montecarlopropagation(circ, wop, 0.5)) <: Tuple{typeof(wop),Bool}
    wop = wrapcoefficients(PauliString(nq, :Z, rand(1:nq)))
    @test typeof(montecarlopropagation(circ, wop, ones(nparams) * π)) <: Tuple{typeof(wop),Bool}

    # TODO Tests including noise channel
end