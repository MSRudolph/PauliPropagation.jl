### plotting.jl
##
# Plotting functionality for visualizing Pauli evolution trees.
# Supports GraphViz output and potentially other visualization backends.
##
###

# Import necessary functions
import ..PauliPropagation: propagate, PauliString, PauliSum, wrapcoefficients, unwrapcoefficients
using JSON
"""
    export_to_graphviz(filename::String="pauli_evolution_tree.dot"; 
                      show_coefficients::Bool=true,
                      node_shape::String="ellipse",
                      edge_style::String="solid")

Export the current evolution tree to a GraphViz DOT file.
"""
function export_to_graphviz(filename::String="pauli_evolution_tree.dot";
    show_coefficients::Bool=true,
    node_shape::String="ellipse",
    edge_style::String="solid")

    open(filename, "w") do io
        println(io, "digraph PauliEvolutionTree {")
        println(io, "    rankdir=TB;")  # Top to Bottom layout
        println(io, "    node [shape=$node_shape, fontname=\"Arial\"];")
        println(io, "    edge [fontname=\"Arial\", fontsize=10];")
        println(io, "")

        # Write nodes
        println(io, "    // Nodes")
        for (node_id, node) in EVOLUTION_TREE
            label = node.pauli_string
            if !isnothing(node.gate_applied)
                label *= "\\n($(node.gate_applied))"
            end
            println(io, "    \"$node_id\" [label=\"$label\", fontsize=12];")
        end

        println(io, "")
        println(io, "    // Edges")
        # Write edges
        for edge in EVOLUTION_EDGES
            if show_coefficients
                label = edge.coefficient
                println(io, "    \"$(edge.parent_id)\" -> \"$(edge.child_id)\" [label=\"$label\", style=$edge_style];")
            else
                println(io, "    \"$(edge.parent_id)\" -> \"$(edge.child_id)\" [style=$edge_style];")
            end
        end

        println(io, "}")
    end

    println("Evolution tree exported to $filename")
    println("To visualize, run: dot -Tpng $filename -o pauli_tree.png")
    return filename
end

"""
    export_to_json(filename::String="pauli_evolution_tree.json")

Export the current evolution tree to JSON format for other visualization tools.
"""
function export_to_json(filename::String="pauli_evolution_tree.json")
    try
        # Try to load JSON package
        JSON = Base.require(@__MODULE__, :JSON)
    catch
        # Fallback if JSON package is not available
        println("Warning: JSON package not available. Install it using: Pkg.add(\"JSON\")")
        println("Exporting as basic text format instead...")

        open(filename, "w") do io
            println(io, "# Pauli Evolution Tree Data")
            println(io, "# Nodes:")
            for (node_id, node) in EVOLUTION_TREE
                println(io, "Node: $node_id, Pauli: $(node.pauli_string), Gate: $(node.gate_applied)")
            end
            println(io, "\n# Edges:")
            for edge in EVOLUTION_EDGES
                println(io, "Edge: $(edge.parent_id) -> $(edge.child_id), Coeff: $(edge.coefficient), Gate: $(edge.gate)")
            end
        end

        println("Evolution tree exported as text to $filename")
        return filename
    end

    # Convert tree data to JSON-friendly format
    tree_data = Dict(
        "nodes" => [
            Dict(
                "id" => node.id,
                "pauli_string" => node.pauli_string,
                "gate_applied" => node.gate_applied
            ) for node in values(EVOLUTION_TREE)
        ],
        "edges" => [
            Dict(
                "parent_id" => edge.parent_id,
                "child_id" => edge.child_id,
                "coefficient" => edge.coefficient,
                "gate" => edge.gate
            ) for edge in EVOLUTION_EDGES
        ]
    )

    open(filename, "w") do io
        JSON.print(io, tree_data, 4)  # Use indent=4 as per user rules
    end

    println("Evolution tree exported to $filename")
    return filename
end

"""
    print_tree_summary()

Print a summary of the current evolution tree.
"""
function print_tree_summary()
    println("=== Pauli Evolution Tree Summary ===")
    println("Number of nodes: $(length(EVOLUTION_TREE))")
    println("Number of edges: $(length(EVOLUTION_EDGES))")

    # Count gates
    gate_counts = Dict{String,Int}()
    for edge in EVOLUTION_EDGES
        gate_counts[edge.gate] = get(gate_counts, edge.gate, 0) + 1
    end

    println("\nGates applied:")
    for (gate, count) in gate_counts
        println("  $gate: $count times")
    end

    # Show root nodes (nodes with no parents)
    root_nodes = [node_id for node_id in keys(EVOLUTION_TREE)
                  if !any(edge.child_id == node_id for edge in EVOLUTION_EDGES)]

    println("\nRoot nodes: $(length(root_nodes))")
    for root_id in root_nodes
        root_node = EVOLUTION_TREE[root_id]
        println("  $root_id: $(root_node.pauli_string)")
    end

    # Show leaf nodes (nodes with no children)
    leaf_nodes = [node_id for node_id in keys(EVOLUTION_TREE)
                  if !any(edge.parent_id == node_id for edge in EVOLUTION_EDGES)]

    println("\nLeaf nodes: $(length(leaf_nodes))")
    for leaf_id in leaf_nodes
        leaf_node = EVOLUTION_TREE[leaf_id]
        println("  $leaf_id: $(leaf_node.pauli_string)")
    end
end

"""
    visualize_tree(output_format::String="graphviz", filename::String="")

Main function for visualizing the evolution tree.
Supports "graphviz", "json", and "summary" formats.
"""
function visualize_tree(output_format::String="graphviz", filename::String="")
    if isempty(EVOLUTION_TREE)
        println("Warning: Evolution tree is empty. Run propagation with PauliTreeTracker coefficients first.")
        return
    end

    if output_format == "graphviz"
        default_filename = "pauli_evolution_tree.dot"
        export_to_graphviz(isempty(filename) ? default_filename : filename)
    elseif output_format == "json"
        default_filename = "pauli_evolution_tree.json"
        export_to_json(isempty(filename) ? default_filename : filename)
    elseif output_format == "summary"
        print_tree_summary()
    else
        throw(ArgumentError("Unsupported output format: $output_format. Use 'graphviz', 'json', or 'summary'."))
    end
end

"""
    propagate_with_tree_tracking(circ, pstr::PauliString, thetas=nothing; kwargs...)

Convenience function that wraps coefficients in PauliTreeTracker and runs propagation.
Returns both the result and exports the tree visualization.
"""
function propagate_with_tree_tracking(circ, pstr::PauliString, thetas=nothing;
    export_format::String="graphviz",
    export_filename::String="",
    reset_tree_first::Bool=true,
    kwargs...)

    if reset_tree_first
        reset_tree!()
    end

    # Wrap the coefficient in PauliTreeTracker
    tracked_pstr = wrapcoefficients(pstr, PauliTreeTracker)

    # Add the initial node
    pstr_str = format_pauli_string(pstr)

    for (pstr_key, coeff) in PauliSum(tracked_pstr)
        add_node!(coeff.node_id, pstr_str, nothing)
        break  # Only one term for PauliString
    end

    # Run propagation
    result = propagate(circ, tracked_pstr, thetas; kwargs...)

    # Export visualization
    visualize_tree(export_format, export_filename)

    # Return unwrapped result
    return unwrapcoefficients(result)
end

"""
    propagate_with_tree_tracking(circ, psum::PauliSum, thetas=nothing; kwargs...)

Convenience function for PauliSum that wraps coefficients in PauliTreeTracker and runs propagation.
"""
function propagate_with_tree_tracking(circ, psum::PauliSum, thetas=nothing;
    export_format::String="graphviz",
    export_filename::String="",
    reset_tree_first::Bool=true,
    kwargs...)

    if reset_tree_first
        reset_tree!()
    end

    # Wrap the coefficients in PauliTreeTracker
    tracked_psum = wrapcoefficients(psum, PauliTreeTracker)

    # Add the initial nodes
    for (pstr, coeff) in tracked_psum
        pstr_str = inttostring(pstr, psum.nqubits)
        add_node!(coeff.node_id, pstr_str, nothing)
    end

    # Run propagation
    result = propagate(circ, tracked_psum, thetas; kwargs...)

    # Export visualization
    visualize_tree(export_format, export_filename)

    # Return unwrapped result
    return unwrapcoefficients(result)
end