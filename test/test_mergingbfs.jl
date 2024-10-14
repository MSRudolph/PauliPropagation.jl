using Random


function numericalPP(nq, nl, W, min_abs_coeff)

    symbs = [:I for _ in 1:nq]
    symbs[round(Int, nq / 2)] = :Z   # as symbol. Also works but is slower.

    obsint = symboltoint(symbs)  # for performance we work with bitoperations

    topo = bricklayertopology(nq; periodic=false)
    circ = hardwareefficientcircuit(nq, nl; topology=topo)
    fastcirc = tofastgates(circ)
    m = countparameters(fastcirc)

    Random.seed!(42)
    thetas = randn(m)

    dnum = mergingbfs(fastcirc, obsint, thetas; max_weight=W, min_abs_coeff=min_abs_coeff)

    return evalagainstzero(dnum) # expectation
end



function hybridPP(nq, nl, W, min_abs_coeff, max_freq)

    symbs = [:I for _ in 1:nq]
    symbs[round(Int, nq / 2)] = :Z   # as symbol. Also works but is slower.

    obsint = symboltoint(symbs)  # for performance we work with bitoperations


    topo = bricklayertopology(nq; periodic=false)
    circ = hardwareefficientcircuit(nq, nl; topology=topo)
    fastcirc = tofastgates(circ)
    m = countparameters(fastcirc)

    Random.seed!(42)
    thetas = randn(m)

    dhyb = mergingbfs(fastcirc, obsint, NumericPathProperties(1.0), thetas; max_weight=W, max_freq=max_freq, min_abs_coeff=min_abs_coeff)

    return evalagainstzero(dhyb)
end


function surrogatePP(nq, nl, W, max_freq)

    symbs = [:I for _ in 1:nq]
    symbs[round(Int, nq / 2)] = :Z   # as symbol. Also works but is slower.

    obsint = symboltoint(symbs)  # for performance we work with bitoperations

    topo = bricklayertopology(nq; periodic=false)
    circ = hardwareefficientcircuit(nq, nl; topology=topo)
    fastcirc = tofastgates(circ)
    m = countparameters(fastcirc)

    Random.seed!(42)
    thetas = randn(m)

    dsym = mergingbfs(circ, operatortopathdict(obsint), zeros(m); max_weight=W, max_freq=max_freq)

    final_nodes = collect(pth.coeff for (obs, pth) in zerofilter(dsym))
    final_eval_node = PauliGateNode(parents=final_nodes, trig_inds=zeros(Int, length(final_nodes)), signs=ones(length(final_nodes)), param_idx=1, cummulative_value=0.0)
    resetnodes(final_eval_node)
    resetnodes(final_eval_node)
    eval_list = gettraceevalorder(final_eval_node, zeros(m))

    return expectation(eval_list, thetas)
end
