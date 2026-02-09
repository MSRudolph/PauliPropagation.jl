# Test for PathProperties.jl
using Test

@testset "Test PathProperties" begin

    cases = [(8, 5, 1.), (11, 4, 0.1)]

    for (i, (nq, j, c)) in enumerate(cases)

        pstr = PauliString(nq, :X, j, PauliFreqTracker(c))

        # Test the PathProiperties coefficients
        @test pstr.coeff.coeff == c

        # Test multiplication
        a = 2
        new_path = pstr.coeff * a
        @test new_path.coeff == a * c
        
    end
end    
