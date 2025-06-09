# Visualization Example
# This script demonstrates how to use the PauliTreeTracker for visualizing
# Pauli string evolution during propagation

using PauliPropagation

# Create a simple circuit with Pauli rotation gates
function create_example_circuit()
    # Create a simple circuit: RX gate followed by RZ gate
    circ = [
        PauliRotation(:X, 1),
        PauliRotation(:Y, 1),
        PauliRotation(:Z, 1),
    ]
    return circ
end

# Create a two-qubit circuit example
function create_two_qubit_circuit()
    # Create a two-qubit circuit with various rotations
    circ = [
        PauliRotation(:Z, 1),    # RZ on qubit 1
        PauliRotation(:Z, 2),    # RZ on qubit 2  
    ]
    return circ
end

# Create a circuit with both Clifford and Pauli rotation gates
function create_mixed_circuit()
    # Create a circuit mixing Clifford and Pauli rotation gates
    circ = [
        PauliRotation(:X, 1),    # RX rotation (branches)
        CliffordGate(:H, 1),     # Hadamard gate (no branching)
        PauliRotation(:Z, 1),    # RZ rotation (branches)
        CliffordGate(:Z, 1),     # Z gate (no branching)
    ]
    return circ
end

function run_one_qubit_example()
    println("=== One-Qubit Example ===")

    # Create a simple 1-qubit Pauli string: Z_1
    nqubits = 1
    pstr = PauliString(nqubits, :Z, 1, 1.0)
    println("Initial Pauli string: $pstr")

    # Create example circuit
    circ = create_example_circuit()
    println("Circuit: RX(θ₁) -> RY(θ₂) -> RZ(θ₃)")

    # Set some parameter values
    thetas = [π / 4, π / 6, π / 8]  # θ₁ = π/4, θ₂ = π/6, θ₃ = π/8

    # Reset tree storage before starting
    reset_tree!()

    # Run propagation with tree tracking
    println("\nRunning propagation with tree tracking...")
    result = propagate_with_tree_tracking(
        circ, pstr, thetas;
        export_format="summary",
        reset_tree_first=true
    )

    println("\nFinal result:")
    println(result)

    # Export to GraphViz format
    println("\nExporting to GraphViz...")
    visualize_tree("graphviz", "pauli_evolution_1qubit.dot")

    # Export to JSON format
    println("Exporting to JSON...")
    visualize_tree("json", "pauli_evolution_1qubit.json")

    # use the normal propagation to verify the result
    println("\nRunning propagation without tree tracking...")
    result = propagate(circ, pstr, thetas)
    println("Final result:")
    println(result)
end

function run_two_qubit_example()
    println("\n=== Two-Qubit Example ===")

    # Create a 2-qubit PauliSum with both XX and ZZ strings
    nqubits = 2

    # Create PauliSum and add terms using the correct method
    psum = PauliSum(nqubits)
    add!(psum, [:X, :X], [1, 2], 1.0)  # XX string with coefficient 1.0
    add!(psum, [:Z, :Z], [1, 2], 0.5)  # ZZ string with coefficient 0.5

    println("Initial Pauli sum: XX + 0.5*ZZ")
    println(psum)

    # Create two-qubit circuit
    circ = create_two_qubit_circuit()
    println("Circuit: RZ₁(θ₁) -> RZ₂(θ₂)")

    # Set some parameter values
    thetas = [π / 3, π / 4]  # Different angles for each gate

    # Reset tree storage before starting
    reset_tree!()

    # Run propagation with tree tracking
    println("\nRunning propagation with tree tracking...")
    result = propagate_with_tree_tracking(
        circ, psum, thetas;
        export_format="summary",
        reset_tree_first=true
    )

    println("\nFinal result:")
    println(result)

    # Export to GraphViz format
    println("\nExporting to GraphViz...")
    visualize_tree("graphviz", "pauli_evolution_2qubit.dot")

    # Export to JSON format
    println("Exporting to JSON...")
    visualize_tree("json", "pauli_evolution_2qubit.json")

    println("\nTwo-qubit example completed!")
    println("To visualize the tree, run:")
    println("  dot -Tpng pauli_evolution_2qubit.dot -o pauli_tree_2qubit.png")

    # use the normal propagation to verify the result
    println("\nRunning propagation without tree tracking...")
    result = propagate(circ, psum, thetas)
    println("Final result:")
    println(result)
end

function run_mixed_circuit_example()
    println("\n=== Mixed Circuit Example (Clifford + Pauli Rotations) ===")

    # Create a simple 1-qubit Pauli string: X_1
    nqubits = 1
    pstr = PauliString(nqubits, :X, 1, 1.0)
    println("Initial Pauli string: $pstr")

    # Create mixed circuit
    circ = create_mixed_circuit()
    println("Circuit: RX(θ₁) -> H -> RZ(θ₂) -> Z")

    # Set some parameter values (only for the parametrized gates)
    thetas = [π / 4, π / 6]  # θ₁ = π/4, θ₂ = π/6

    # Reset tree storage before starting
    reset_tree!()

    # Run propagation with tree tracking
    println("\nRunning propagation with tree tracking...")
    result = propagate_with_tree_tracking(
        circ, pstr, thetas;
        export_format="summary",
        reset_tree_first=true
    )

    println("\nFinal result:")
    println(result)

    # Export to GraphViz format
    println("\nExporting to GraphViz...")
    visualize_tree("graphviz", "pauli_evolution_mixed.dot")

    # Export to JSON format
    println("Exporting to JSON...")
    visualize_tree("json", "pauli_evolution_mixed.json")

    println("\nMixed circuit example completed!")
    println("To visualize the tree, run:")
    println("  dot -Tpng pauli_evolution_mixed.dot -o pauli_tree_mixed.png")

    # use the normal propagation to verify the result
    println("\nRunning propagation without tree tracking...")
    result = propagate(circ, pstr, thetas)
    println("Final result:")
    println(result)
end

function main()
    println("=== Pauli Evolution Visualization Examples ===")

    # Run one-qubit example
    # run_one_qubit_example()

    # Run two-qubit example
    # run_two_qubit_example()

    # Run mixed circuit example
    run_mixed_circuit_example()

    println("\n=== All Examples Completed ===")
    println("Generated files:")
    println("  - pauli_evolution_1qubit.dot/json (1-qubit Z -> RX,RY,RZ)")
    println("  - pauli_evolution_2qubit.dot/json (2-qubit XX -> RY₁,RX₂,RZ₁,RY₂)")
    println("  - pauli_evolution_mixed.dot/json (1-qubit X -> RX,H,RZ,Z)")
end

# Run the example
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end