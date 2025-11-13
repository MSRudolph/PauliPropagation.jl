using Test

@testset "pauliprod Tests" begin
    """Test the product of PauliStrings and PauliSums."""

    # X * Y = iZ
    nq = 2
    qind = 1
    pstr1 = PauliString(nq, :X, qind, 1)
    pstr2 = PauliString(nq, :Y, qind, 1.5)
    pstr3 = pauliprod(pstr1, pstr2)
    @test pstr3 == PauliString(nq, :Z, qind, pstr1.coeff * pstr2.coeff * 1im)
    @test pauliprod(pstr2, pstr1).coeff == -1 * pstr3.coeff

    # Z * I = Z
    nq = 5
    qind = 4
    pstr1 = PauliString(nq, :Z, qind, 1.2 + 0.2im)
    pstr2 = PauliString(nq, :I, qind, 1im)
    pstr3 = pauliprod(pstr1, pstr2)
    @test pstr3 == PauliString(nq, :Z, qind, pstr1.coeff * pstr2.coeff)
    @test pauliprod(pstr2, pstr1).coeff == pstr3.coeff

    # Z * I = Z
    nq = 7
    qind = 7
    pstr1 = PauliString(nq, :I, qind, -1im)
    pstr2 = PauliString(nq, :X, qind, 0)
    pstr3 = pauliprod(pstr1, pstr2)
    @test pstr3 == PauliString(nq, :X, qind, pstr1.coeff * pstr2.coeff)
    @test pauliprod(pstr2, pstr1).coeff == pstr3.coeff


    # X * X = I
    nq = 13
    qind = 2
    pstr1 = PauliString(nq, :X, qind, 1im + 0.5)
    pstr2 = PauliString(nq, :X, qind, -0.2)
    pstr3 = pauliprod(pstr1, pstr2)
    @test pstr3 == PauliString(nq, :I, qind, pstr1.coeff * pstr2.coeff)
    @test pauliprod(pstr2, pstr1).coeff == pstr3.coeff


    # XY * ZY = -iYI
    nq = 17
    qinds = [3, 14]
    pstr1 = PauliString(nq, [:X, :Y], qinds, 0.3)
    pstr2 = PauliString(nq, [:Z, :Y], qinds, 5)
    pstr3 = pauliprod(pstr1, pstr2)
    @test pstr3 == PauliString(nq, [:Y, :I], qinds, pstr1.coeff * pstr2.coeff * -1im * 1)
    @test pauliprod(pstr2, pstr1).coeff == -1 * pstr3.coeff


    # I * I = I
    nq = 32
    qind = 16
    pstr1 = PauliString(nq, :I, qind, 2.1)
    pstr2 = PauliString(nq, :I, qind, 3im + 0.1)
    pstr3 = pauliprod(pstr1, pstr2)
    @test pstr3 == PauliString(nq, :I, qind, pstr1.coeff * pstr2.coeff)
    @test pauliprod(pstr2, pstr1).coeff == pstr3.coeff

    # XY * ZY = -iYI
    nq = 33
    qinds = [1, 29]
    pstr1 = PauliString(nq, [:X, :I], qinds, 2im)
    pstr2 = PauliString(nq, [:I, :Y], qinds, -3)
    pstr3 = pauliprod(pstr1, pstr2)
    @test pstr3 == PauliString(nq, [:X, :Y], qinds, pstr1.coeff * pstr2.coeff)
    @test pauliprod(pstr2, pstr1).coeff == pstr3.coeff

    # ZX * XZ = YY
    nq = 65
    qinds = [62, 31]
    pstr1 = PauliString(nq, [:Z, :X], qinds, 1im)
    pstr2 = PauliString(nq, [:X, :Z], qinds, π)
    pstr3 = pauliprod(pstr1, pstr2)
    @test pstr3 == PauliString(nq, [:Y, :Y], qinds, pstr1.coeff * pstr2.coeff)
    @test pauliprod(pstr2, pstr1).coeff == pstr3.coeff
end

@testset "LinearAlgebra" begin
    for nq in 1:10
        λ = rand()
        @test tr(λ * PauliString(nq, :I, rand(1:nq))) == λ * 2.0^nq
        @test tr(λ * PauliString(nq, :X, rand(1:nq))) == 0.0
        @test tr(λ * PauliString(nq, :Y, rand(1:nq))) == 0.0
        @test tr(λ * PauliString(nq, :Z, rand(1:nq))) == 0.0

        @test trace(λ * PauliString(nq, :I, rand(1:nq))) == λ * 2.0^nq
        @test trace(λ * PauliString(nq, :X, rand(1:nq))) == 0.0
        @test trace(λ * PauliString(nq, :Y, rand(1:nq))) == 0.0
        @test trace(λ * PauliString(nq, :Z, rand(1:nq))) == 0.0        
    end

    nq = 42
    letters = [rand([:X, :Y, :Z]) for _ in 1:nq]
    pstr = PauliString(nq, letters, 1:nq)
    @test tr(pstr) == 0.0
    @test trace(pstr) == 0.0

    psum = pstr + PauliString(nq, :I, 1) / 137.
    @test tr(psum) == 2.0^nq / 137.
    @test trace(psum) == 2.0^nq / 137.
end

# TODO: Tests for PauliSum
# TODO: Test for commutator