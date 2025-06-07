# Pauli Evolution Visualization

This submodule provides tree-like visualization capabilities for tracking how Pauli strings evolve during propagation through quantum circuits.

## Overview

The visualization system tracks the genealogy of Pauli strings as they transform through gates, creating a tree structure that shows:
- **Nodes**: Individual Pauli strings (e.g., "XXZ", "IYI")  
- **Edges**: Transformation coefficients (e.g., "cos(θ)", "-sin(θ)")

## Key Components

### `PauliTreeTracker`
A specialized `PathProperties` subtype that tracks:
- Numerical coefficient (`coeff`)
- Unique node ID (`node_id`) 
- Parent-child relationships (`parent_id`)

### Tree Storage
Global storage containers:
- `EVOLUTION_TREE`: Dictionary mapping node IDs to `TreeNode` objects
- `EVOLUTION_EDGES`: Vector of `TreeEdge` objects representing transformations

### Visualization Functions
- `visualize_tree()`: Main visualization function
- `export_to_graphviz()`: Exports to GraphViz DOT format
- `export_to_json()`: Exports to JSON format
- `print_tree_summary()`: Prints tree statistics

## Usage Example

```julia
using PauliPropagation

# Create a circuit and initial Pauli string
circ = [PauliRotation(:X, 1), PauliRotation(:Z, 1)]
pstr = PauliString(1, :Z, 1, 1.0)
thetas = [π/4, π/6]

# Run propagation with tree tracking
result = propagate_with_tree_tracking(circ, pstr, thetas)

# Visualize the evolution tree
visualize_tree("graphviz", "my_tree.dot")  # GraphViz format
visualize_tree("json", "my_tree.json")     # JSON format  
visualize_tree("summary")                  # Print summary

# Generate PNG image (requires GraphViz installed)
# shell> dot -Tpng my_tree.dot -o my_tree.png
```

## Tree Structure

For a Pauli rotation gate that doesn't commute with a Pauli string, the transformation creates two branches:

```
Original Pauli (e.g., "Z")
├── cos(θ) → Same Pauli ("Z") 
└── ±sin(θ) → New Pauli (e.g., "Y")
```

## Output Formats

### GraphViz (.dot)
- Nodes show Pauli strings and the gate that created them
- Edges show the multiplication coefficients
- Can be rendered to PNG, SVG, PDF using GraphViz tools

### JSON
- Machine-readable format for custom visualization tools
- Contains complete node and edge information
- Formatted with 4-space indentation

### Summary
- Text-based overview of tree statistics
- Shows number of nodes, edges, gates applied
- Lists root and leaf nodes

## Advanced Usage

### Manual Tree Construction
```julia
# Reset tree before starting
reset_tree!()

# Manually add nodes and edges
add_node!("node1", "XYZ", "RX(θ)")
add_edge!("parent1", "node1", "cos(θ)", "RX(θ)")

# Use with existing PauliSum containing PauliTreeTracker coefficients
tracked_psum = wrapcoefficients(my_psum, PauliTreeTracker)
result = propagate(circ, tracked_psum, thetas)
```

### Customization
- Node shapes and edge styles can be customized in GraphViz export
- Tree traversal functions can be extended for custom analysis
- Merge behavior can be modified by overriding the `+` operator

## Limitations

- Currently focuses on demonstration rather than scalability
- Tree tracking adds computational overhead
- Memory usage grows with circuit depth and Pauli string count
- Merge and truncation events are simplified in the current implementation

## Dependencies

- **UUIDs**: For generating unique node identifiers
- **JSON** (optional): For JSON export functionality. Install with `Pkg.add("JSON")`
- **GraphViz** (external): For rendering DOT files to images 