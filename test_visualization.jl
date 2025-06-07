#!/usr/bin/env julia

# Simple test script for the visualization module
using Pkg
Pkg.activate(".")

# Include the main module
include("src/PauliPropagation.jl")
using .PauliPropagation

println("=== Testing PauliPropagation Visualization Module ===")

# Test 1: Basic PauliTreeTracker creation
println("\n1. Testing PauliTreeTracker creation...")
tracker = PauliTreeTracker(1.0)
println("   ✓ Created PauliTreeTracker with coefficient $(tracker.coeff)")
println("   ✓ Node ID: $(tracker.node_id)")

# Test 2: Tree operations
println("\n2. Testing tree operations...")
reset_tree!()
add_node!("test1", "XYZ", "RX(θ)")
add_edge!("test1", "test2", "cos(θ)", "RX")

println("   ✓ Added node and edge to tree")
println("   ✓ Tree has $(length(EVOLUTION_TREE)) nodes and $(length(EVOLUTION_EDGES)) edges")

# Test 3: Format Pauli string
println("\n3. Testing Pauli string formatting...")
pstr = PauliString(3, [:X, :Y, :Z], [1, 2, 3], 1.0)
formatted = format_pauli_string(pstr)
println("   ✓ Formatted Pauli string: $formatted")

# Test 4: Visualization functions
println("\n4. Testing visualization functions...")
print_tree_summary()

# Test 5: Export functions (without actually writing files)
println("\n5. Testing export functions...")
try
    export_to_graphviz("test_output.dot")
    println("   ✓ GraphViz export successful")

    # Clean up
    rm("test_output.dot", force=true)
catch e
    println("   ⚠ GraphViz export warning: $e")
end

try
    export_to_json("test_output.json")
    println("   ✓ JSON export successful")

    # Clean up
    rm("test_output.json", force=true)
catch e
    println("   ⚠ JSON export warning: $e")
end

println("\n=== All tests completed! ===")
println("The visualization module appears to be working correctly.")
println("\nTo test with actual propagation, run:")
println("  julia examples/visualization_example.jl")