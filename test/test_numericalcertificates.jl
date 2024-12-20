@testset "numerical certificate Tests" begin
    nq = rand(1:100)

    nl = rand(1:10)
    topo = bricklayertopology(nq; periodic=false)

    circ = efficientsu2circuit(nq, nl; topology=topo)
    nparams = countparameters(circ)
    ngates = length(circ)


    # Numerical Coefficient
    pstr = PauliString(nq, :X, nq)
    @test isa(estimateaverageerror(circ, pstr, 100), Float64)
    pstr = PauliString(nq, :Y, nq)
    @test isa(estimateaverageerror(circ, pstr, 100, 0.87236582), Float64)

    # Weight Truncation
    pstr = PauliString(nq, :Z, nq)
    @test estimateaverageerror(circ, pstr, 100, π; max_weight=0) >= estimateaverageerror(circ, pstr, 100, π; max_weight=nq) == 0.0
    @test estimateaverageerror(circ, pstr, 100, ones(nparams) * 1.23; max_weight=0) >= estimateaverageerror(circ, pstr, 100, ones(nparams) * 1.23; max_weight=nq) == 0.0

    # Frequency Truncation
    wrapped_pstr = wrapcoefficients(PauliString(nq, :Z, nq), PauliFreqTracker)
    @test estimateaverageerror(circ, wrapped_pstr, 100, π; max_freq=0) >= estimateaverageerror(circ, pstr, 100, π; max_freq=nq) == 0.0
    @test estimateaverageerror(circ, wrapped_pstr, 100, ones(nparams) * 1.23; max_freq=0) >= estimateaverageerror(circ, pstr, 100, ones(nparams) * 1.23; max_freq=nq) == 0.0

    # Small-angle Truncation
    wrapped_pstr = wrapcoefficients(PauliString(nq, :Z, nq), PauliFreqTracker)
    @test estimateaverageerror(circ, wrapped_pstr, 100, π; max_sins=0) >= estimateaverageerror(circ, pstr, 100, π; max_sins=nq) == 0.0
    @test estimateaverageerror(circ, wrapped_pstr, 100, ones(nparams) * 1.23; max_sins=0) >= estimateaverageerror(circ, pstr, 100, ones(nparams) * 1.23; max_sins=nq) == 0.0

end