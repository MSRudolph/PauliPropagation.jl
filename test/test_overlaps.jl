# Test File for stateoverlap.jl
# 
# This file contains unit tests for the following functionalities:
# - Edge cases for overlapwithones.
# - TODO(YT): add tests for other functions in stateoverlap.jl

using Test

test_cases = [
    "twoZs" => (4, (2, 4), 1.0),
    "oneZ" => (5, (2), -1.0),
    "noZ" => (4, (3), 1.0),
    "noones" => (4, (), 1.0),
]

@testset "Parameterized Tests stateoverlaps" begin

    for (name, (nq, indices, expected)) in test_cases
        @testset "$name" begin

            # SubTest the case where the PauliString is the identity
            pstr = PauliString(nq, :I, nq)
            @test overlapwithones(pstr, indices) == 1.0

        end

        @testset "$name" begin

            # SubTest the case where the PauliString is a Z operator
            pstr = PauliString(nq, [:Z, :Z], [2, 4])
            @test overlapwithones(pstr, indices) == expected

        end

        @testset "$name" begin

            # SubTest the case where the PauliString contains an X/Y operator
            pstr = PauliString(nq, [:X, :Z, :Z], [1, 2, 4])
            @test overlapwithones(pstr, indices) == 0.

            pstr = PauliString(nq, [:Y], [4])
            @test overlapwithones(pstr, indices) == 0.

        end
    end

end