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
    CliffordGate(symbol::Symbol, qinds::Union{AbstractArray, Tuple, Base.Generator})

Constructor for a `CliffordGate` acting on the qubits `qinds`. 
Converts the types of `qinds` to the correct types for `CliffordGate`.
"""
function CliffordGate(symbol::Symbol, qinds::Union{AbstractArray,Tuple,Base.Generator})
    return CliffordGate(symbol, collect(qinds))
end

# TODO: verify that these are all correct
const _default_clifford_map = Dict(
    :H => [(1, 0x00), (1, 0x03), (-1, 0x02), (1, 0x01)],
    :X => [(1, 0x00), (1, 0x01), (-1, 0x02), (-1, 0x03)],
    :Y => [(1, 0x00), (-1, 0x01), (1, 0x02), (1, 0x03)],
    :Z => [(1, 0x00), (-1, 0x01), (-1, 0x02), (1, 0x03)],
    :S => [(1, 0x00), (-1, 0x02), (1, 0x01), (1, 0x03)],
    :CNOT => [(1, 0x00), (1, 0x05), (1, 0x06), (1, 0x03), (1, 0x04), (1, 0x01), (1, 0x02), (1, 0x07), (1, 0x0b), (1, 0x0e), (-1, 0x0d), (1, 0x08), (1, 0x0f), (-1, 0x0a), (1, 0x09), (1, 0x0c)],
    :ZZpihalf => [(1, 0x00), (1, 0x0e), (-1, 0x0d), (1, 0x03), (1, 0x0b), (1, 0x05), (1, 0x06), (1, 0x08), (-1, 0x07), (1, 0x09), (1, 0x0a), (-1, 0x04), (1, 0x0c), (1, 0x02), (-1, 0x01), (1, 0x0f)],
    :SWAP => [
        (1, 0x00), (1, 0x04), (1, 0x08), (1, 0x0c), (1, 0x01), (1, 0x05),
        (1, 0x09), (1, 0x0d), (1, 0x02), (1, 0x06), (1, 0x0a), (1, 0x0e),
        (1, 0x03), (1, 0x07), (1, 0x0b), (1, 0x0f)
    ],
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
    createcliffordmap(gate_relations::Dict)

Create a Clifford gate map from a dictionary of gate relations which can then be pushed to the global `clifford_map`.
`gate_relations` is a dictionary with pairs like `(:X, :X) => (-1, :Z, :X)`,
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
        Tuple{Int,typeof(gate_keys[1]).parameters...}
    }(undef, length(gate_keys))
    for (i, idx) in enumerate(order_indices)
        reordered_gate_vals[idx+1] = gate_relations[gate_keys[i]]
    end

    mapped_gate = Vector{Tuple{Int,UInt8}}(undef, length(gate_keys))
    for (i, v) in enumerate(reordered_gate_vals)
        mapped_gate[i] = v[1], symboltoint(collect(v[2:end]))
    end

    return mapped_gate
end

## Apply CliffordGate helpers
"""
    applywithmap(gate::CliffordGate, pstr::PauliStringType, coefficient, map_array)

Apply a `CliffordGate` to an integer Pauli string and a coefficient 
using the a `map_array` corresponding to the `CliffordGate`.
"""
function applywithmap(gate::CliffordGate, pstr::PauliStringType, coefficient, map_array; kwargs...)
    qinds = gate.qinds

    lookup_int = _extractlookupop(pstr, qinds)
    sign, partial_pstr = map_array[lookup_int+1]  # +1 because Julia is 1-indexed and lookup_int is 0-indexed
    pstr = _insertnewpaulis!(pstr, partial_pstr, qinds)

    coefficient = _multiplysign(coefficient, sign)
    return pstr, coefficient
end

function _extractlookupop(lookup_int::PauliStringType, qinds)
    partial_pstr = typeof(lookup_int)(0)
    for ii in eachindex(qinds)
        partial_pstr = setpauli(partial_pstr, getpauli(lookup_int, qinds[ii]), ii)
    end
    return partial_pstr
end

function _insertnewpaulis!(pstr::PauliStringType, partial_pstr::PauliStringType, qinds)
    for ii in eachindex(qinds)
        pstr = setpauli(pstr, getpauli(partial_pstr, ii), qinds[ii])
    end
    return pstr
end

function _multiplysign(coefficient, sign)
    return coefficient * sign
end