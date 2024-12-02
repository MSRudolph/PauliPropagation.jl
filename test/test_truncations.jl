using Random
using Test

function brickcircuit(seed)
    nq = 8
    nl = 4
    op = PauliString(nq, :Z, round(Int, nq / 2))

    topo = bricklayertopology(nq; periodic=false)
    circ = hardwareefficientcircuit(nq, nl; topology=topo)

    Random.seed!(seed)
    thetas = randn(length(circ))

    return circ, op, thetas
end

@testset "Truncate damping coefficients Tests" begin
    """Test the truncations by damped coefficients."""
    seed = 42
    circ, op, thetas = brickcircuit(seed)

    W = Inf
    min_abs_coeff = 0.0
    evolved_p = propagate(
        circ, op, thetas;
        max_weight=W, min_abs_coeff=min_abs_coeff
    )
    expected_expval = overlapwithzero(evolved_p)

    gamma = 0.0
    truncategamma = (pstr, coeff) -> truncatedampingcoeff(
        pstr, coeff, gamma, min_abs_coeff
    )
    evolved_p = propagate(
        circ, op, thetas;
        max_weight=W, min_abs_coeff=min_abs_coeff,
        customtruncatefn=truncategamma
    )
    # \gamma=0 == zero dissipation
    @test isapprox(overlapwithzero(evolved_p), expected_expval)

    gamma = 0.01
    min_abs_coeff = 1e-5
    truncategamma = (pstr, coeff) -> truncatedampingcoeff(
        pstr, coeff, gamma, min_abs_coeff
    )
    evolved_p = propagate(
        circ, op, thetas;
        max_weight=W, min_abs_coeff=min_abs_coeff,
        customtruncatefn=truncategamma
    )
    # \gamma=0.1 \approx zero dissipation
    #TODO: is there another way to test dissipation?
    @test isapprox(overlapwithzero(evolved_p), expected_expval, rtol=1e-3)

end