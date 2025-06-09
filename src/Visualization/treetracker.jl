### treetracker.jl
##
# PauliTreeTracker type and methods for tracking genealogy during propagation.
# It records parent-child relationships and edge coefficients during gate application.
##
###

using UUIDs

# Import necessary functions from parent modules
import ..PauliPropagation: PathProperties, PauliSum, PauliString, PauliStringType, PauliRotation, MaskedPauliRotation, CliffordGate
import ..PauliPropagation: _tomaskedpaulirotation, paulitype, getnewpaulistring, commutes, set!, add!
import ..PauliPropagation: inttostring, symboltoint, getpauli, setpauli, splitapply, applytoall!, apply
import ..PauliPropagation: _applysin, _applycos

"""
    TreeNode

Represents a node in the Pauli evolution tree.
Contains the Pauli string representation and a unique identifier.
"""
struct TreeNode
    id::String
    pauli_string::String
    gate_applied::Union{String,Nothing}  # The gate that created this node (if any)
end

"""
    TreeEdge

Represents an edge in the Pauli evolution tree.
Contains the coefficient that was applied during the transformation.
"""
struct TreeEdge
    parent_id::String
    child_id::String
    coefficient::String  # String representation of the coefficient (e.g., "cos(θ₁)", "-sin(θ₁)")
    gate::String
end

"""
    PauliTreeTracker(coeff::Number, node_id::String, parent_id::Union{String, Nothing})

Wrapper type for numerical coefficients in Pauli propagation that tracks the genealogy tree.
Each PauliTreeTracker corresponds to a node in the evolution tree.
"""
mutable struct PauliTreeTracker{T<:Number} <: PathProperties
    coeff::T
    node_id::String
    parent_id::Union{String,Nothing}

    # Constructor
    function PauliTreeTracker(coeff::T, node_id::String, parent_id::Union{String,Nothing}=nothing) where {T<:Number}
        new{T}(coeff, node_id, parent_id)
    end
end

"""
    PauliTreeTracker(coeff::Number)

Constructor for `PauliTreeTracker` from only a coefficient.
Creates a new unique node ID and no parent.
"""
function PauliTreeTracker(coeff::Number)
    node_id = string(uuid4())[1:8]  # Use first 8 characters of UUID for readability
    return PauliTreeTracker(float(coeff), node_id, nothing)
end

# Global storage for the evolution tree
const EVOLUTION_TREE = Dict{String,TreeNode}()
const EVOLUTION_EDGES = Vector{TreeEdge}()

"""
    reset_tree!()

Reset the global evolution tree storage.
"""
function reset_tree!()
    empty!(EVOLUTION_TREE)
    empty!(EVOLUTION_EDGES)
end

"""
    add_node!(node_id::String, pauli_string::String, gate_applied::Union{String, Nothing}=nothing)

Add a node to the evolution tree.
"""
function add_node!(node_id::String, pauli_string::String, gate_applied::Union{String,Nothing}=nothing)
    println("Adding node: $node_id, $pauli_string, $gate_applied")
    EVOLUTION_TREE[node_id] = TreeNode(node_id, pauli_string, gate_applied)
end

"""
    add_edge!(parent_id::String, child_id::String, coefficient::String, gate::String)

Add an edge to the evolution tree.
"""
function add_edge!(parent_id::String, child_id::String, coefficient::String, gate::String)
    push!(EVOLUTION_EDGES, TreeEdge(parent_id, child_id, coefficient, gate))
end

"""
    create_child_tracker(parent::PauliTreeTracker, coeff_expression::String, gate_name::String)

Create a new child tracker from a parent during gate application.
"""
function create_child_tracker(parent::PauliTreeTracker, coeff_value::Number, coeff_expression::String, gate_name::String)
    child_id = string(uuid4())[1:8]
    child_tracker = PauliTreeTracker(coeff_value, child_id, parent.node_id)

    # Add edge to the tree
    add_edge!(parent.node_id, child_id, coeff_expression, gate_name)

    return child_tracker
end

"""
    format_pauli_string(pstr)

Convert a Pauli string to a readable string format like "XXZ".
"""
function format_pauli_string(pstr::PauliStringType, nqubits::Int)
    return inttostring(pstr, nqubits)
end

function format_pauli_string(pstr::PauliString)
    return inttostring(pstr.term, pstr.nqubits)
end

function format_pauli_string(pstr)
    # Fallback for other types
    return string(pstr)
end

# Override the addition operation to handle tree merging
function Base.:+(path1::PauliTreeTracker, path2::PauliTreeTracker)
    # When merging, we create a new node that represents the sum
    # In practice, you might want to keep track of which paths were merged
    merged_coeff = path1.coeff + path2.coeff
    merged_id = string(uuid4())[1:8]

    # Create merged tracker - for simplicity, we'll use path1's parent
    merged_tracker = PauliTreeTracker(merged_coeff, merged_id, path1.parent_id)

    # Add edges showing the merge (optional for visualization)
    if !isnothing(path1.parent_id) && !isnothing(path2.parent_id)
        add_edge!(path1.node_id, merged_id, "merge_1", "MERGE")
        add_edge!(path2.node_id, merged_id, "merge_2", "MERGE")
    end

    return merged_tracker
end

### Specialized methods for gate applications

"""
    splitapply(gate::MaskedPauliRotation, pstr::PauliStringType, coeff::PauliTreeTracker, theta; nqubits::Int, kwargs...)

Specialized splitapply for PauliTreeTracker that tracks the tree evolution.
This creates two child nodes with cos and sin coefficients.
"""
function splitapply(gate::MaskedPauliRotation, pstr::PauliStringType, coeff::PauliTreeTracker, theta; nqubits::Int, kwargs...)
    # Get the gate name for labeling - extract first symbol from gate.symbols
    gate_symbol = isempty(gate.symbols) ? "?" : string(gate.symbols[1])
    gate_name = "R$(gate_symbol)"

    # Add the current node to the tree if not already there
    pstr_str = format_pauli_string(pstr, nqubits)


    # Create cos coefficient child
    cos_coeff_value = coeff.coeff * cos(theta)
    cos_multiplier = cos(theta)
    cos_child = create_child_tracker(coeff, cos_coeff_value, string(round(cos_multiplier, digits=3)), gate_name)
    add_node!(cos_child.node_id, pstr_str, gate_name)

    # Get new Pauli string and sign for sin coefficient
    new_pstr, sign = getnewpaulistring(gate, pstr)
    sin_coeff_value = coeff.coeff * sin(theta) * sign
    sin_multiplier = sin(theta) * sign
    sin_child = create_child_tracker(coeff, sin_coeff_value, string(round(sin_multiplier, digits=3)), gate_name)

    # Add the new Pauli string node
    new_pstr_str = format_pauli_string(new_pstr, nqubits)
    add_node!(sin_child.node_id, new_pstr_str, gate_name)

    return pstr, cos_child, new_pstr, sin_child
end

"""
    _applysin(pth::PauliTreeTracker, theta, sign=1; kwargs...)

Apply sin coefficient while tracking the genealogy.
"""
function _applysin(pth::PauliTreeTracker, theta, sign=1; kwargs...)
    sin_coeff = pth.coeff * sin(theta) * sign
    return PauliTreeTracker(sin_coeff, pth.node_id, pth.parent_id)
end

"""
    _applycos(pth::PauliTreeTracker, theta, sign=1; kwargs...)

Apply cos coefficient while tracking the genealogy.
"""
function _applycos(pth::PauliTreeTracker, theta, sign=1; kwargs...)
    cos_coeff = pth.coeff * cos(theta) * sign
    return PauliTreeTracker(cos_coeff, pth.node_id, pth.parent_id)
end

"""
    applytoall!(gate::PauliRotation, theta, psum::PauliSum{TT,PauliTreeTracker{T}}, aux_psum; kwargs...)

Specialized applytoall! for PauliRotation gates with PauliSum containing PauliTreeTracker coefficients.
This tracks every gate application in the evolution tree.
"""
function applytoall!(gate::PauliRotation, theta, psum::PauliSum{TT,PauliTreeTracker{T}}, aux_psum; kwargs...) where {TT<:PauliStringType,T<:Number}
    # Convert PauliRotation to MaskedPauliRotation for efficiency
    gate = _tomaskedpaulirotation(gate, paulitype(psum))

    # Loop over all Pauli strings and their coefficients in the Pauli sum
    for (pstr, coeff) in psum
        pstr_str = format_pauli_string(pstr, psum.nqubits)
        println("pstr_str: $pstr_str")
        if commutes(gate, pstr)
            # If the gate commutes, create a new child node and edge
            pstr_str = format_pauli_string(pstr, psum.nqubits)
            gate_symbol = isempty(gate.symbols) ? "?" : string(gate.symbols[1])
            gate_name = "R$(gate_symbol)"

            # Create new child tracker (coefficient stays the same for commuting gates)
            new_child = create_child_tracker(coeff, coeff.coeff, "1", gate_name)
            add_node!(new_child.node_id, pstr_str, gate_name)

            # Update the coefficient in the sum with the new child tracker
            set!(psum, pstr, new_child)
            continue
        end

        # Apply the gate and track the split
        pstr, coeff1, new_pstr, coeff2 = splitapply(gate, pstr, coeff, theta; nqubits=psum.nqubits, kwargs...)

        # Set the coefficient of the original Pauli string
        set!(psum, pstr, coeff1)

        # Set the coefficient of the new Pauli string in the aux_psum
        set!(aux_psum, new_pstr, coeff2)
    end

    return
end

"""
    applytoall!(gate::CliffordGate, theta, psum::PauliSum{TT,PauliTreeTracker{T}}, aux_psum; kwargs...)

Specialized applytoall! for CliffordGate with PauliSum containing PauliTreeTracker coefficients.
Clifford gates deterministically transform Pauli strings without branching, so we create a single child node for each transformation.
"""
function applytoall!(gate::CliffordGate, theta, psum::PauliSum{TT,PauliTreeTracker{T}}, aux_psum; kwargs...) where {TT<:PauliStringType,T<:Number}
    # Loop over all Pauli strings and their coefficients in the Pauli sum
    for (pstr, coeff) in psum
        # Apply the Clifford gate to get the new Pauli string and coefficient
        new_pstr, new_coeff_value = apply(gate, pstr, coeff.coeff; kwargs...)

        # Format the gate name for display
        gate_name = string(gate.symbol)

        # Create a new child tracker for the transformed Pauli string
        edge_num = new_coeff_value / coeff.coeff
        new_child = create_child_tracker(coeff, new_coeff_value, string(round(edge_num, digits=3)), gate_name)

        # Add the new node to the tree
        new_pstr_str = format_pauli_string(new_pstr, psum.nqubits)
        add_node!(new_child.node_id, new_pstr_str, gate_name)

        # Set the new Pauli string and its tracker in aux_psum
        # (Note: Clifford gates create non-overlapping Pauli strings so we can use set!)
        set!(aux_psum, new_pstr, new_child)
    end

    # Empty the original psum since everything was moved to aux_psum
    empty!(psum)

    return
end