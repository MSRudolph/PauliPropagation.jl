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
