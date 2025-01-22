using Random


function numericalPP(nq, nl, W, min_abs_coeff)

    pstr = PauliString(nq, :Z, round(Int, nq / 2))

    topo = bricklayertopology(nq; periodic=false)
    circ = hardwareefficientcircuit(nq, nl; topology=topo)

    m = countparameters(circ)

    Random.seed!(42)
    thetas = randn(m)

    dnum = propagate(circ, pstr, thetas; max_weight=W, min_abs_coeff=min_abs_coeff)

    return overlapwithzero(dnum) # expectation
end



function hybridPP(nq, nl, W, min_abs_coeff, max_freq)

    pstr = PauliString(nq, :Z, round(Int, nq / 2))

    wrapped_pstr = wrapcoefficients(pstr, PauliFreqTracker)

    topo = bricklayertopology(nq; periodic=false)
    circ = hardwareefficientcircuit(nq, nl; topology=topo)

    m = countparameters(circ)

    Random.seed!(42)
    thetas = randn(m)

    dhyb = propagate(circ, wrapped_pstr, thetas; max_weight=W, max_freq=max_freq, min_abs_coeff=min_abs_coeff)

    return overlapwithzero(dhyb)
end


function surrogatePP(nq, nl, W, max_freq)

    pstr = PauliString(nq, :Z, round(Int, nq / 2))

    wrapped_pstr = wrapcoefficients(pstr, NodePathProperties)

    topo = bricklayertopology(nq; periodic=false)
    circ = hardwareefficientcircuit(nq, nl; topology=topo)

    m = countparameters(circ)

    Random.seed!(42)
    thetas = randn(m)

    dsym = propagate(circ, wrapped_pstr; max_weight=W, max_freq=max_freq)
    zerofilter!(dsym)  # Filter the nodes that you find relevant
    evaluate!(dsym, thetas)

    return overlapwithzero(dsym)
end

@testset begin
    nq = 6
    nl = 3

    topology = rectangletopology(2, 3; periodic=true)
    circ = efficientsu2circuit(nq, nl; topology=topology)

    nparams = countparameters(circ)

    for max_weight in (1, 3, 6)
        thetas = randn(nparams)

        pstr = PauliString(nq, rand([:I, :X, :Y, :Z]), rand(1:nq))
        dnum = propagate(circ, pstr, thetas; min_abs_coeff=0, max_weight=max_weight)

        wrapped_pstr = wrapcoefficients(pstr, PauliFreqTracker)
        dhyb = propagate(circ, wrapped_pstr, thetas; min_abs_coeff=0, max_weight=max_weight)

        surrogate_pstr = wrapcoefficients(pstr, NodePathProperties)
        dsym = propagate(circ, surrogate_pstr; max_weight=max_weight)
        evaluate!(dsym, thetas)

        @test overlapwithzero(dnum) ≈ overlapwithzero(dhyb) ≈ overlapwithzero(dsym)
        @test overlapwithplus(dnum) ≈ overlapwithplus(dhyb) ≈ overlapwithplus(dsym)
    end


    nl = 5
    circ = efficientsu2circuit(nq, nl; topology=topology)
    nparams = countparameters(circ)
    thetas = randn(nparams)

    pstr = PauliString(nq, rand([:X, :Y, :Z]), rand(1:nq))
    dnum = propagate(circ, pstr, thetas; min_abs_coeff=1e-3)

    wrapped_pstr = wrapcoefficients(pstr, PauliFreqTracker)
    dhyb = propagate(circ, wrapped_pstr, thetas; min_abs_coeff=1e-3)

    @test overlapwithzero(dnum) ≈ overlapwithzero(dhyb)
    @test overlapwithplus(dnum) ≈ overlapwithplus(dhyb)

end
