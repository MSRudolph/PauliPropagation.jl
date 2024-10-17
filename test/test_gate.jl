using Test

@testset "create_gate_map Tests" begin
    """Test gate map function from a user-defined gate."""
    # CNOT
    CNOT_relations = Dict(  # from IBM
        (:I, :I) => (1, :I, :I),
        (:I, :X) => (1, :I, :X),
        (:I, :Y) => (1, :Z, :Y),
        (:I, :Z) => (1, :Z, :Z),
        (:X, :I) => (1, :X, :X),
        (:X, :X) => (1, :X, :I),
        (:X, :Y) => (1, :Y, :Z),
        (:X, :Z) => (-1, :Y, :Y),
        (:Y, :I) => (1, :Y, :X),
        (:Y, :X) => (1, :Y, :I),
        (:Y, :Y) => (-1, :X, :Z),
        (:Y, :Z) => (1, :X, :Y),
        (:Z, :I) => (1, :Z, :I),
        (:Z, :X) => (1, :Z, :X),
        (:Z, :Y) => (1, :I, :Y),
        (:Z, :Z) => (1, :I, :Z),
    )    
    mapped_CNOT = create_clifford_map(CNOT_relations)
    # Check we get the same result 
    @test mapped_CNOT == default_clifford_map[:CNOT]

    # H
    H_relations = Dict(  # from IBM
        (:I,) => (1, :I),
        (:X,) => (1, :Z),
        (:Y,) => (-1, :Y),
        (:Z,) => (1, :X),
    )
    mapped_H = create_clifford_map(H_relations)
    # Check we get the same result
    @test mapped_H == default_clifford_map[:H]
end
