## Test the utils for PauliAlgebra/datatypes.jl module.

using Test

const symbols = [:I, :X, :Y, :Z]

@testset "Int to Symbol" begin

    # Test the conversion from integer to symbol
    for (i, symbol) in enumerate(symbols)
        @test inttosymbol(i-1) == symbol # i-1 because integers are 0-indexed
    end

    # Test the bounds error for converting integers outside [0, 3]
    @test_throws BoundsError inttosymbol(4)
    @test_throws BoundsError inttosymbol(-1)

    # Test multiqubit conversion. Note that the order of the symbols is reversed
    # compared to the bitstring to int conversion.
    @test inttosymbol(1, 2) == [:X, :I]
    @test inttosymbol(4, 2) == [:I, :X]
end

@testset "Symbol to Int" begin
    
    # Test the conversion from a symbol to integer
    for i in 0:3
        @test symboltoint(inttosymbol(i)) == i
    end

    # Test 2-qubit conversion
    nqubits = 2
    for i in 4:15
        @test symboltoint(inttosymbol(i, nqubits)) == i
    end

end

@testset "Get Pauli" begin

    # For single qubit pauli, the integer representation is the same.
    for i in 0:3
        @test getpauli(i, 1) == i
    end

    pstr = PauliString(3, [:X], [2])
    # Test for `pstr` of `PauliString` type.
    @test getpauli(pstr, 1) == 0
    @test getpauli(pstr, 2) == 1

    # Test for `pstr` of `PauliStringType` type.
    @test getpauli(pstr.term, 1) == 0
    @test getpauli(pstr.term, 2) == 1
end

@testset "Get Paulis" begin
    
    # Test with `PauliString` type
    pstr = PauliString(3, [:X, :Y, :Z], [1, 2, 3])
    @test getpaulis(pstr, [1, 2, 3]) == pstr
    
    @test getpaulis(pstr, [1]) == PauliString(1, [:X], [1])
    @test getpaulis(pstr, [1, 2]) == PauliString(2, [:X, :Y], [1, 2])

    # Test consistency with integer representation
    rand_int = 13
    @test getpaulis(rand_int, [1, 2]) == 13

    # Test with `PauliStringType`
    nqubits = 4
    pstr = PauliString(nqubits, :X, 2)
    inds = [1, 3]
    @test getpaulis(pstr.term, inds) == 0
    inds = [3, 1]
    @test getpaulis(pstr.term, inds) == 0
    inds = [1, 2]
    @test getpaulis(pstr.term, inds) == 4
    # if the indices are not sorted, `getpaulis` will reutrn a pauli
    # integer representation according to the unsorted indices
    inds = [2, 1]
    @test getpaulis(pstr.term, inds) == 1

    # Test with inttosymbol
    inds = [2, 4]
    symbols = [:X, :Z]
    pstr = PauliString(nqubits, symbols, inds)

    @test inttosymbol(
        getpaulis(pstr.term, inds), length(inds)
    ) == symbols

    # flip `qinds` will return reversed symbols
    @test inttosymbol(
        getpaulis(pstr.term, [4, 2]), length(inds)
    ) == [symbols[2], symbols[1]]    

end

@testset "Set Pauli" begin
    nqubits = 4
    # Test pstr with `PauliString` type
    pstr = PauliString(nqubits, [:X, :Y, :Z], [1, 2, 3])
    @test setpauli(pstr, :Z, 1) == PauliString(nqubits, [:Z, :Y, :Z], [1, 2, 3])

    # Test pstr with `PauliStringType`
    pstr = PauliString(nqubits, :X, 2)
    expected_pstr = PauliString(nqubits, [:Z, :X], [1, 2])
    @test setpauli(pstr.term, :Z, 1) == expected_pstr.term

    # Test target_pauli as integer
    @test setpauli(pstr, 3, 1) == expected_pstr
    @test setpauli(pstr.term, 3, 1) == expected_pstr.term

end

@testset "Set Paulis" begin
    
    nqubits = 4
    pstr = PauliString(nqubits, [:X, :Z], [1, 4])

    # Test pstr with `PauliString` type and target_pauli as `Integer`
    transformed_pstr = setpaulis(pstr, 2, [1, 3])
    expected_pstr = PauliString(4, [:Y, :Z], [1, 4])
    @test transformed_pstr == expected_pstr

    # Test pstr with `PauliString` type and target_pauli as `Symbol`
    transformed_pstr = setpaulis(pstr, [:Y, :X], [1, 3])
    expected_pstr = PauliString(4, [:Y, :X, :Z], [1, 3, 4])    
    @test transformed_pstr == expected_pstr    

    # Test pstr with `PauliStringType` and target_pauli as `Integer`
    transformed_pstr = setpaulis(pstr.term, symboltoint([:Y, :X]), [1, 3])
    @test transformed_pstr == expected_pstr.term

    expected_pstr = PauliString(4, [:Y, :X, :Z], [3, 1, 4])  
    transformed_pstr = setpaulis(pstr.term, [:Y, :X], [3, 1])
    @test transformed_pstr == expected_pstr.term    
end

@testset "Int to PauliString" begin
    
    # Test the conversion from integer to PauliString
    @test inttopstr(0, 2) == PauliString(2, [:I, :I], [1, 2])
    @test inttopstr(1, 2) == PauliString(2, [:X, :I], [1, 2])
    @test inttopstr(2, 2) == PauliString(2, [:Y, :I], [1, 2])
    @test inttopstr(12, 2) == PauliString(2, [:I, :Z], [1, 2])

    # Test consistency with integer representation
    rand_int = 13
    @test inttopstr(rand_int, 2).term == rand_int

end


