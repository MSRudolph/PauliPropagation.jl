using Test

@testset "Test create Clifford gates" begin
    """Test gate map function from a user-defined gate."""
    # CNOT
    CNOT_relations = Dict(
        (:I, :I) => (:I, :I, 1),
        (:I, :X) => (:I, :X, 1),
        (:I, :Y) => (:Z, :Y, 1),
        (:I, :Z) => (:Z, :Z, 1),
        (:X, :I) => (:X, :X, 1),
        (:X, :X) => (:X, :I, 1),
        (:X, :Y) => (:Y, :Z, 1),
        (:X, :Z) => (:Y, :Y, -1),
        (:Y, :I) => (:Y, :X, 1),
        (:Y, :X) => (:Y, :I, 1),
        (:Y, :Y) => (:X, :Z, -1),
        (:Y, :Z) => (:X, :Y, 1),
        (:Z, :I) => (:Z, :I, 1),
        (:Z, :X) => (:Z, :X, 1),
        (:Z, :Y) => (:I, :Y, 1),
        (:Z, :Z) => (:I, :Z, 1),
    )
    mapped_CNOT = createcliffordmap(CNOT_relations)
    @test mapped_CNOT == clifford_map[:CNOT]

    clifford_map[:CNOT2] = mapped_CNOT
    reset_clifford_map!()
    @test clifford_map == PauliPropagation._default_clifford_map

    # H
    H_relations = Dict(
        (:I,) => (:I, 1),
        (:X,) => (:Z, 1),
        (:Y,) => (:Y, -1),
        (:Z,) => (:X, 1),
    )
    mapped_H = createcliffordmap(H_relations)
    @test mapped_H == clifford_map[:H]

end


@testset "Test transposing Clifford maps" begin
    for (gate, map) in clifford_map
        transposed_map = transposecliffordmap(map)
        @test length(transposed_map) == length(map)
        @test transposecliffordmap(transposed_map) == map
    end
end


@testset "Test concatenating Clifford maps" begin
    circuit = [CliffordGate(:H, [2]), CliffordGate(:CNOT, [1, 2]), CliffordGate(:H, [2])]
    @test concatenatecliffordmaps(circuit) == clifford_map[:CZ]

    circuit = [CliffordGate(:H, [2]), CliffordGate(:CZ, [1, 2]), CliffordGate(:H, [2])]
    @test concatenatecliffordmaps(circuit) == clifford_map[:CNOT]

    circuit = [CliffordGate(:H, [1]), CliffordGate(:Z, [1]), CliffordGate(:H, [1])]
    @test concatenatecliffordmaps(circuit) == clifford_map[:X]

    circuit = [CliffordGate(:SX, [1]), CliffordGate(:SX, [1])]
    @test concatenatecliffordmaps(circuit) == clifford_map[:X]

    circuit = [CliffordGate(:S, [1]), CliffordGate(:S, [1])]
    @test concatenatecliffordmaps(circuit) == clifford_map[:Z]

end