using PauliPropagation
using BenchmarkTools
using Random


function timingnumericalPP()
    nq = 8
    nl = 4
    W = Inf
    min_abs_coeff = 0

    op = PauliString(nq, :Z, round(Int, nq / 2))
    opsum = PauliSum(nq, op)

    topo = bricklayertopology(nq; periodic=false)
    # topo = get2dtopology(4, 4)
    circ = hardwareefficientcircuit(nq, nl; topology=topo)
    fastcirc = tofastgates(circ)
    m = length(fastcirc)

    Random.seed!(42)
    thetas = randn(m)

    res1 = mergingbfs(fastcirc, op, thetas; max_weight=W, min_abs_coeff=min_abs_coeff)
    res2 = mergingbfs(fastcirc, opsum, thetas; max_weight=W, min_abs_coeff=min_abs_coeff)
    @show overlapwithzero(res1), overlapwithzero(res2)
    @btime mergingbfs($fastcirc, $op, $thetas; max_weight=$W, min_abs_coeff=$min_abs_coeff)
    @btime mergingbfs($fastcirc, $opsum, $thetas; max_weight=$W, min_abs_coeff=$min_abs_coeff)

    return
end

# 55.345 ms (354 allocations: 3.79 MiB)
timingnumericalPP()