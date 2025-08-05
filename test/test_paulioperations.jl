using Test

@testset "pauliprod Tests" begin
    """Test the product of PauliStrings and PauliSums."""

    # Test with PauliStrings
    @testset "PauliString Product" begin
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

    @testset "PauliSums" begin

        psum1 = PauliSum(3)
        add!(psum1, [:Z, :X], [1, 2])

        psum2 = PauliSum(3)
        add!(psum2, [:Y, :Y], [1, 3])
        add!(psum2, [:X, :Y], [1, 2])

        # Psum * I = Psum
        @test pauliprod(psum2, PauliSum(3, [:I], [1])) == psum2

        # Expected product
        psum_expected = PauliSum(3)
        add!(psum_expected, [:Y, :Z], [1, 2], -1.)
        add!(psum_expected, [:X, :X, :Y], [1, 2, 3], -1im)
        @test pauliprod(psum1, psum2) == psum_expected

        # Test real coefficients
        psum2 = PauliSum(3, [:Y, :Y], [1, 3], 2.0)
        psum_expected = PauliSum(3)
        add!(psum_expected, [:Y, :Z], [1, 2], -2.0)
        @test pauliprod(psum1, psum2) == psum_expected
    end
end
