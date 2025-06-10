# Visualization Example
# This script demonstrates how to use the PauliTreeTracker for visualizing
# Pauli string evolution during propagation

using PauliPropagation
using Random

function create_three_qubit_circuit()
    # create a hardware efficient ansatz circuit
    nq = 3
    nl = 2
    topology = bricklayertopology(nq; periodic=false)
    circuit = hardwareefficientcircuit(nq, nl; topology=topology)
    nparams = countparameters(circuit)
    return circuit, nparams
end
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
        CliffordGate(:X, 1),
        PauliRotation(:X, 1),    # RX rotation (branches)
        CliffordGate(:H, 1),     # Hadamard gate (no branching)
        PauliRotation(:Z, 1),    # RZ rotation (branches)
        CliffordGate(:Z, 1),     # Z gate (no branching)
    ]
    return circ
end

# Create a noisy circuit with Pauli noise channels
function create_pauli_noise_circuit()
    # Create a circuit with various Pauli noise channels
    circ = [
        PauliRotation(:X, 1),      # RX rotation
        DepolarizingNoise(1),      # Depolarizing noise on qubit 1
        PauliRotation(:Y, 2),      # RY rotation on qubit 2
        PauliXNoise(2),            # Pauli-X noise on qubit 2
        PauliRotation(:Z, 1),      # RZ rotation on qubit 1
        PauliZNoise(1),            # Pauli-Z (dephasing) noise on qubit 1
    ]
    return circ
end

# Create a circuit with amplitude damping noise
function create_amplitude_damping_circuit()
    # Create a circuit with amplitude damping noise
    circ = [
        PauliRotation(:X, 1),      # RX rotation on qubit 1
        AmplitudeDampingNoise(1),  # Amplitude damping on qubit 1
        PauliRotation(:Y, 2),      # RY rotation on qubit 2
        AmplitudeDampingNoise(2),  # Amplitude damping on qubit 2
        PauliRotation(:Z, 1),      # RZ rotation on qubit 1
        AmplitudeDampingNoise(1),  # Another amplitude damping on qubit 1
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
    psum = PauliSum(nqubits)
    add!(psum, [:X], [1], 1.0)
    add!(psum, [:Y], [1], 1.0)
    add!(psum, [:Z], [1], 1.0)

    println("Initial Pauli sum: $psum")

    # Create mixed circuit
    circ = create_mixed_circuit()
    println("Circuit: X -> RX(θ₁) -> H -> RZ(θ₂) -> Z")

    # Set some parameter values (only for the parametrized gates)
    thetas = [π / 4, π / 6]  # θ₁ = π/4, θ₂ = π/6

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
    visualize_tree("graphviz", "pauli_evolution_mixed.dot")

    # Export to JSON format
    println("Exporting to JSON...")
    visualize_tree("json", "pauli_evolution_mixed.json")

    println("\nMixed circuit example completed!")
    println("To visualize the tree, run:")
    println("  dot -Tpng pauli_evolution_mixed.dot -o pauli_tree_mixed.png")

    # use the normal propagation to verify the result
    println("\nRunning propagation without tree tracking...")
    result = propagate(circ, psum, thetas)
    println("Final result reference:")
    println(result)
end

function run_mixed_multiple_qubit_example()
    println("\n=== Mixed Multiple Qubit Example (Clifford + Pauli Rotations) ===")
    nqubits = 3
    psum = PauliSum(nqubits)
    add!(psum, [:I, :I, :X], [1, 2, 3], 1.0)

    println("Initial Pauli sum: $psum")

    # Create mixed circuit
    circ, nparams = create_three_qubit_circuit()

    # Set some parameter values (only for the parametrized gates)
    Random.seed!(42)
    thetas = randn(nparams) * 0.5

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
    visualize_tree("graphviz", "pauli_evolution_mixed_multiple_qubit.dot")
    # run normal propagation to verify the result
    println("\nRunning propagation without tree tracking...")
    result = propagate(circ, psum, thetas)
    println("Final result reference:")
    println(result)
end

function run_pauli_noise_example()
    println("\n=== Pauli Noise Circuit Example ===")

    # Create a 2-qubit PauliSum with X and Z strings
    nqubits = 2
    psum = PauliSum(nqubits)
    add!(psum, [:X], [1], 1.0)    # X string on qubit 1
    add!(psum, [:Z], [2], 0.5)    # Z string on qubit 2

    println("Initial Pauli sum: X₁ + 0.5*Z₂")
    println(psum)

    # Create Pauli noise circuit
    circ = create_pauli_noise_circuit()
    println("Circuit: RX₁(θ₁) -> DepolarizingNoise₁(p₁) -> RY₂(θ₂) -> PauliXNoise₂(p₂) -> RZ₁(θ₃) -> PauliZNoise₁(p₃)")

    # Set parameter values: [rotation angles, noise strengths]
    # First 3 are rotation angles, last 3 are noise strengths
    thetas = [π / 4, π / 6, π / 8, 0.1, 0.05, 0.15]

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
    visualize_tree("graphviz", "pauli_noise_evolution.dot")

    # Export to JSON format
    println("Exporting to JSON...")
    visualize_tree("json", "pauli_noise_evolution.json")

    println("\nPauli noise example completed!")
    println("To visualize the tree, run:")
    println("  dot -Tpng pauli_noise_evolution.dot -o pauli_noise_tree.png")

    # Run normal propagation to verify the result
    println("\nRunning propagation without tree tracking...")
    result_verify = propagate(circ, psum, thetas)
    println("Final result reference:")
    println(result_verify)
end

function run_amplitude_damping_example()
    println("\n=== Amplitude Damping Circuit Example ===")

    # Create a 2-qubit PauliSum with various Pauli strings
    nqubits = 2
    psum = PauliSum(nqubits)
    add!(psum, [:X], [1], 1.0)    # X string on qubit 1 (will be damped)
    add!(psum, [:Y], [2], 0.8)    # Y string on qubit 2 (will be damped)
    add!(psum, [:Z], [1], 0.6)    # Z string on qubit 1 (will split)

    println("Initial Pauli sum: X₁ + 0.8*Y₂ + 0.6*Z₁")
    println(psum)

    # Create amplitude damping circuit
    circ = create_amplitude_damping_circuit()
    println("Circuit: RX₁(θ₁) -> AmplitudeDamping₁(γ₁) -> RY₂(θ₂) -> AmplitudeDamping₂(γ₂) -> RZ₁(θ₃) -> AmplitudeDamping₁(γ₃)")

    # Set parameter values: [rotation angles, damping rates]
    # First 3 are rotation angles, last 3 are damping rates (γ values)
    thetas = [π / 3, π / 4, π / 6, 0.2, 0.1, 0.15]

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
    visualize_tree("graphviz", "amplitude_damping_evolution.dot")

    # Export to JSON format
    println("Exporting to JSON...")
    visualize_tree("json", "amplitude_damping_evolution.json")

    println("\nAmplitude damping example completed!")
    println("To visualize the tree, run:")
    println("  dot -Tpng amplitude_damping_evolution.dot -o amplitude_damping_tree.png")

    # Run normal propagation to verify the result
    println("\nRunning propagation without tree tracking...")
    result_verify = propagate(circ, psum, thetas)
    println("Final result reference:")
    println(result_verify)
end

function main()
    println("=== Pauli Evolution Visualization Examples ===")

    # Run one-qubit example
    # run_one_qubit_example()

    # Run two-qubit example
    # run_two_qubit_example()

    # Run mixed circuit example
    # run_mixed_circuit_example()

    # Run mixed multiple qubit example
    # run_mixed_multiple_qubit_example()

    # Run Pauli noise example
    run_pauli_noise_example()

    # Run amplitude damping example
    run_amplitude_damping_example()

    println("\n=== All Examples Completed ===")
    println("Generated files:")
    println("  - pauli_evolution_1qubit.dot/json (1-qubit Z -> RX,RY,RZ)")
    println("  - pauli_evolution_2qubit.dot/json (2-qubit XX -> RY₁,RX₂,RZ₁,RY₂)")
    println("  - pauli_evolution_mixed.dot/json (1-qubit X -> RX,H,RZ,Z)")
    println("  - pauli_noise_evolution.dot/json (2-qubit noisy circuit with Pauli noise)")
    println("  - amplitude_damping_evolution.dot/json (2-qubit circuit with amplitude damping)")
end

# Run the example
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end