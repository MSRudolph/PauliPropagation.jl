@testset "numerical certificate Tests" begin
    nq = rand(1:100)

    nl = 4
    topo = bricklayertopology(nq; periodic=false)

    circ = efficientsu2circuit(nq, nl; topology=topo)
    nparams = countparameters(circ)
    ngates = length(circ)


    # Numerical Coefficient
    pstr = PauliString(nq, :X, rand(1:nq))
    @test isa(estimateaverageerror(circ, pstr, 100), Float64)
    pstr = PauliString(nq, :Y, rand(1:nq))
    @test isa(estimateaverageerror(circ, pstr, 100, 0.87236582), Float64)

    # Weight Truncation
    pstr = PauliString(nq, :Z, rand(1:nq))
    @test estimateaverageerror(circ, pstr, 100, π; max_weight=0) >= estimateaverageerror(circ, pstr, 100, π; max_weight=nq) == 0.0
    @test estimateaverageerror(circ, pstr, 100, ones(nparams) * 1.23; max_weight=0) >= estimateaverageerror(circ, pstr, 100, ones(nparams) * 1.23; max_weight=nq) == 0.0

    # Frequency Truncation
    wop = wrapcoefficients(PauliString(nq, :Z, rand(1:nq)))
    @test estimateaverageerror(circ, wop, 100, π; max_freq=0) >= estimateaverageerror(circ, pstr, 100, π; max_freq=nq) == 0.0
    @test estimateaverageerror(circ, wop, 100, ones(nparams) * 1.23; max_freq=0) >= estimateaverageerror(circ, pstr, 100, ones(nparams) * 1.23; max_freq=nq) == 0.0

    # Small-angle Truncation
    wop = wrapcoefficients(PauliString(nq, :Z, rand(1:nq)))
    @test estimateaverageerror(circ, wop, 100, π; max_sins=0) >= estimateaverageerror(circ, pstr, 100, π; max_sins=nq) == 0.0
    @test estimateaverageerror(circ, wop, 100, ones(nparams) * 1.23; max_sins=0) >= estimateaverageerror(circ, pstr, 100, ones(nparams) * 1.23; max_sins=nq) == 0.0

    # Numerical Coefficient
    pstr = PauliString(nq, :X, rand(1:nq))
    @test typeof(montecarlopropagation(circ, pstr)) <: Tuple{typeof(pstr),Bool}
    pstr = PauliString(nq, :Y, rand(1:nq))
    @test typeof(montecarlopropagation(circ, pstr, 0.5)) <: Tuple{typeof(pstr),Bool}
    pstr = PauliString(nq, :Z, rand(1:nq))
    @test typeof(montecarlopropagation(circ, pstr, ones(nparams) * π)) <: Tuple{typeof(pstr),Bool}

    # NumericPathProperties Coefficient
    wop = wrapcoefficients(PauliString(nq, :X, rand(1:nq)))
    @test typeof(montecarlopropagation(circ, wop)) <: Tuple{typeof(wop),Bool}
    wop = wrapcoefficients(PauliString(nq, :Y, rand(1:nq)))
    @test typeof(montecarlopropagation(circ, wop, 0.5)) <: Tuple{typeof(wop),Bool}
    wop = wrapcoefficients(PauliString(nq, :Z, rand(1:nq)))
    @test typeof(montecarlopropagation(circ, wop, ones(nparams) * π)) <: Tuple{typeof(wop),Bool}

    # TODO Tests including noise channel
end