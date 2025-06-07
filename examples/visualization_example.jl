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

function main()
    println("=== Pauli Evolution Visualization Example ===")

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
    visualize_tree("graphviz", "pauli_evolution_example.dot")

    # Export to JSON format
    println("Exporting to JSON...")
    visualize_tree("json", "pauli_evolution_example.json")

    println("\nExample completed!")
    println("To visualize the tree, run:")
    println("  dot -Tpng pauli_evolution_example.dot -o pauli_tree.png")

    # use the normal propagation to verify the result
    println("\nRunning propagation without tree tracking...")
    result = propagate(circ, pstr, thetas)
    println("Final result:")
    println(result)
end

# Run the example
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end