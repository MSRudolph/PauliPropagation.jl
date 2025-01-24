"""
    CliffordGate(symbol::Symbol, qinds::Vector{Int})

A Clifford gate with the name `symbol` acting on the qubits `qinds`.
`symbol` needs to match any of the implemented Clifford gates in the global `clifford_map`.
"""
struct CliffordGate <: StaticGate
    symbol::Symbol
    qinds::Vector{Int}
end

"""
    CliffordGate(symbol::Symbol, qind::Int)

Constructor for a single-qubit `CliffordGate`.
"""
function CliffordGate(symbol::Symbol, qind::Int)
    return CliffordGate(symbol, [qind])
end

"""
    CliffordGate(symbol::Symbol, qinds::Tuple{Int...})

Constructor for a `CliffordGate` acting on the qubits `qinds`. 
Converts the types of `qinds` to the correct types for `CliffordGate`.
"""
function CliffordGate(symbol::Symbol, qinds::Union{AbstractArray,Tuple,Base.Generator})
    return CliffordGate(symbol, collect(qinds))
end

# TODO: verify that these are all correct
# const _default_clifford_map = Dict(
#     :H => [(1, 0x00), (1, 0x03), (-1, 0x02), (1, 0x01)],
#     :X => [(1, 0x00), (1, 0x01), (-1, 0x02), (-1, 0x03)],
#     :Y => [(1, 0x00), (-1, 0x01), (1, 0x02), (1, 0x03)],
#     :Z => [(1, 0x00), (-1, 0x01), (-1, 0x02), (1, 0x03)],
#     :SX => [(1, 0x00), (1, 0x01), (-1, 0x03), (1, 0x02)],
#     :SY => [(1, 0x00), (1, 0x03), (1, 0x02), (-1, 0x01)],
#     :S => [(1, 0x00), (-1, 0x02), (1, 0x01), (1, 0x03)],
#     :CNOT => [(1, 0x00), (1, 0x05), (1, 0x06), (1, 0x03), (1, 0x04), (1, 0x01), (1, 0x02), (1, 0x07), (1, 0x0b), (1, 0x0e), (-1, 0x0d), (1, 0x08), (1, 0x0f), (-1, 0x0a), (1, 0x09), (1, 0x0c)],
#     :CZ => [(1, 0x00), (1, 0x0d), (1, 0x0e), (1, 0x03), (1, 0x07), (1, 0x0a), (-1, 0x09), (1, 0x04), (1, 0x0b), (-1, 0x06), (1, 0x05), (1, 0x08), (1, 0x0c), (1, 0x01), (1, 0x02), (1, 0x0f)],
#     :ZZpihalf => [(1, 0x00), (1, 0x0e), (-1, 0x0d), (1, 0x03), (1, 0x0b), (1, 0x05), (1, 0x06), (1, 0x08), (-1, 0x07), (1, 0x09), (1, 0x0a), (-1, 0x04), (1, 0x0c), (1, 0x02), (-1, 0x01), (1, 0x0f)],
#     :SWAP => [
#         (1, 0x00), (1, 0x04), (1, 0x08), (1, 0x0c), (1, 0x01), (1, 0x05),
#         (1, 0x09), (1, 0x0d), (1, 0x02), (1, 0x06), (1, 0x0a), (1, 0x0e),
#         (1, 0x03), (1, 0x07), (1, 0x0b), (1, 0x0f)
#     ],
# )
const _default_clifford_map = Dict{Symbol,Vector{Tuple{UInt8,Int64}}}(
    :H => [(0x00, 1), (0x03, 1), (0x02, -1), (0x01, 1)],
    :X => [(0x00, 1), (0x01, 1), (0x02, -1), (0x03, -1)],
    :Y => [(0x00, 1), (0x01, -1), (0x02, 1), (0x03, 1)],
    :Z => [(0x00, 1), (0x01, -1), (0x02, -1), (0x03, 1)],
    :SX => [(0x00, 1), (0x01, 1), (0x03, -1), (0x02, 1)],
    :SY => [(0x00, 1), (0x03, 1), (0x02, 1), (0x01, -1)],
    :S => [(0x00, 1), (0x02, -1), (0x01, 1), (0x03, 1)],
    :CNOT => [(0x00, 1), (0x05, 1), (0x06, 1), (0x03, 1), (0x04, 1), (0x01, 1), (0x02, 1), (0x07, 1), (0x0b, 1), (0x0e, 1), (0x0d, -1), (0x08, 1), (0x0f, 1), (0x0a, -1), (0x09, 1), (0x0c, 1)],
    :CZ => [(0x00, 1), (0x0d, 1), (0x0e, 1), (0x03, 1), (0x07, 1), (0x0a, 1), (0x09, -1), (0x04, 1), (0x0b, 1), (0x06, -1), (0x05, 1), (0x08, 1), (0x0c, 1), (0x01, 1), (0x02, 1), (0x0f, 1)],
    :SWAP => [(0x00, 1), (0x04, 1), (0x08, 1), (0x0c, 1), (0x01, 1), (0x05, 1), (0x09, 1), (0x0d, 1), (0x02, 1), (0x06, 1), (0x0a, 1), (0x0e, 1), (0x03, 1), (0x07, 1), (0x0b, 1), (0x0f, 1)],
    :ZZpihalf => [(0x00, 1), (0x0e, 1), (0x0d, -1), (0x03, 1), (0x0b, 1), (0x05, 1), (0x06, 1), (0x08, 1), (0x07, -1), (0x09, 1), (0x0a, 1), (0x04, -1), (0x0c, 1), (0x02, 1), (0x01, -1), (0x0f, 1)],
)




"""
    clifford_map

Global dictionary of Clifford gates and their action on Pauli strings.
Currently supported Clifford gates are `:H`, `:X`, `:Y`, `:Z`, `:S`, `:CNOT`, `:ZZpihalf`, and `:SWAP`.
If one indexes into the returned arrays with the integer that corresponds to the partial Pauli string,
the returned tuple is `(sign, partial_pstr)` where `sign` is the sign change and `partial_pstr` is the new partial Pauli string.
"""
const clifford_map = deepcopy(_default_clifford_map)

"""
    reset_clifford_map!()

Reset global `clifford_map` to the CLifford gate implemented by default.
"""
function reset_clifford_map!()
    println(
        "Resetting the global clifford_map to the default Clifford gates.\n
        The warning may be ignored."
    )
    global clifford_map = deepcopy(_default_clifford_map)
    return
end

"""
    transposecliffordmap(map_array::Vector{Tuple{UInt8,Int}})

Transpose the Clifford gate `maparray` so that the output map is the inverse of the input map.
For example, `transposecliffordmap(clifford_map[:H])` returns the map for the inverse of the Hadamard gate, which is the same map.
"""
function transposecliffordmap(map_array::Vector{Tuple{UInt8,Int64}})
    original_paulis = [pauli for (pauli, _) in map_array]
    perm_inds = sortperm(original_paulis)
    return Tuple{UInt8,Int}[(UInt8(ii - 1), s) for (ii, (_, s)) in enumerate(map_array)][perm_inds]
end

"""
    createcliffordmap(gate_relations::Dict)

Create a Clifford gate map from a dictionary of gate relations which can then be pushed to the global `clifford_map`.
`gate_relations` is a dictionary with pairs like `(:X, :X) => (:Z, :X, -1)`,
describing the action of the Clifford gate on symbols (including the sign change).
"""
function createcliffordmap(gate_relations::Dict)
    # check that number of qubits are 4 or less (supported by UInt8)
    if maximum([length(k) for k in keys(gate_relations)]) > 4
        throw(ArgumentError(
            "Number of qubits must be 4 or less due to UInt8 type restrictions."
        ))
    end

    gate_keys = collect(keys(gate_relations))
    order_indices = [symboltoint(collect(k)) for k in gate_keys]

    # Initialize arrays for reordered keys and values
    reordered_gate_vals = Vector{
        Tuple{typeof(gate_keys[1]).parameters...,Int}
    }(undef, length(gate_keys))
    for (i, idx) in enumerate(order_indices)
        reordered_gate_vals[idx+1] = gate_relations[gate_keys[i]]
    end

    mapped_gate = Vector{Tuple{UInt8,Int}}(undef, length(gate_keys))
    for (i, v) in enumerate(reordered_gate_vals)
        mapped_gate[i] = symboltoint(collect(v[1:end-1])), v[end]
    end

    return mapped_gate
end


"""
    concatenatecliffordmaps(circuit::Vector{CliffordGate})

Concatenate a circuit of Clifford gates into a single Clifford map.
The length of the map is `4^nq`` where `nq` is the maximum qubit index in the circuit.
The resulting clifford map can be added to the global `clifford_map` with a custom Clifford gate name.
The maximum number of qubits is 4 due to current restrictions of `UInt8`.
Even if all gates only act on one qubit, that qubit index will determine the dimensionality of the map.
"""
function concatenatecliffordmaps(circuit)
    if !(all([isa(gate, CliffordGate) for gate in circuit]))
        throw(ArgumentError("All gates in the circuit must be Clifford gates."))
    end

    max_nq = maximum(maximum(gate.qinds) for gate in circuit)

    if max_nq > 4
        throw(ArgumentError("Number of qubits must be 4 or less due current to UInt8 type restrictions."))
    end

    # max integer to feed into the circuit
    max_integer = 4^max_nq - 1

    psums = [propagate(circuit, PauliString(max_nq, ii, 1.0)) for ii in 0:max_integer]
    if any(length(psum) > 1 for psum in psums)
        throw(ArgumentError("The circuit does not implement a 1 to 1 map of Pauli strings."))
    end
    # pairs of pstr => sign
    pairs = [first(ps) for ps in psums]

    # Convert into a Vector{UInt8, Int} like other clifford maps
    return [(UInt8(pair[1]), Int(pair[2])) for pair in pairs]

end