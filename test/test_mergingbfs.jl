using Random


function numericalPP(nq, nl, W, min_abs_coeff)

    op = PauliString(nq, :Z, round(Int, nq / 2))

    topo = bricklayertopology(nq; periodic=false)
    circ = hardwareefficientcircuit(nq, nl; topology=topo)
    fastcirc = tofastgates(circ)
    m = length(fastcirc)

    Random.seed!(42)
    thetas = randn(m)

    dnum = mergingbfs(fastcirc, op, thetas; max_weight=W, min_abs_coeff=min_abs_coeff)

    return overlapwithzero(dnum) # expectation
end



function hybridPP(nq, nl, W, min_abs_coeff, max_freq)

    op = PauliString(nq, :Z, round(Int, nq / 2))

    d = Dict(op.operator => NumericPathProperties(1.0))  # TODO: Adapt to PauliSum

    topo = bricklayertopology(nq; periodic=false)
    circ = hardwareefficientcircuit(nq, nl; topology=topo)
    fastcirc = tofastgates(circ)
    m = length(fastcirc)

    Random.seed!(42)
    thetas = randn(m)

    dhyb = mergingbfs(fastcirc, d, thetas; max_weight=W, max_freq=max_freq, min_abs_coeff=min_abs_coeff)

    return overlapwithzero(dhyb)
end


function surrogatePP(nq, nl, W, max_freq)

    op = PauliString(nq, :Z, round(Int, nq / 2))

    d = operatortopathdict(op.operator)

    topo = bricklayertopology(nq; periodic=false)
    circ = hardwareefficientcircuit(nq, nl; topology=topo)
    fastcirc = tofastgates(circ)
    m = length(fastcirc)

    Random.seed!(42)
    thetas = randn(m)

    dsym = mergingbfs(circ, d, zeros(m); max_weight=W, max_freq=max_freq)

    final_nodes = collect(pth.coeff for (obs, pth) in zerofilter(dsym))
    final_eval_node = PauliGateNode(parents=final_nodes, trig_inds=zeros(Int, length(final_nodes)), signs=ones(length(final_nodes)), param_idx=1, cummulative_value=0.0)
    resetnodes(final_eval_node)
    resetnodes(final_eval_node)
    eval_list = gettraceevalorder(final_eval_node, zeros(m))

    return expectation(eval_list, thetas)
end
