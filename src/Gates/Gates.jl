"""
Abstract type for gates. 
"""
abstract type Gate end

"""
Abstract type for parametrized gates.
"""
abstract type ParametrizedGate <: Gate end

"""
Abstract type for static gates are not parametrized.
"""
abstract type StaticGate <: Gate end

include("frozengates.jl")
include("pauligates.jl")
include("cliffordgates.jl")
include("noisechannels.jl")


## Interface for transforming gates to potentially optimized gates
"""
    tofastgates(gate::Gate, nqubits::Integer)

Transforms a gate to a potentially faster but more involved gate type when the total number of qubits `nqubits` is known.`
This is currently only for `PauliGate` to `FastPauliGate`.
"""
tofastgates(gate::Gate, nqubits::Integer) = gate


"""
    tofastgates(circ::Vector{G}) where {G<:Gate}

Transforms a circuit in the form of a vector of gates to a vector of potentially faster gates where applicable.
The maximum number of qubits is determined from the gates in the circuit, but they all require a `qinds` field.
"""
function tofastgates(circ::Vector{G}) where {G<:Gate}
    # Find the maximum number of qubits
    nqubits = _getmaxqubits(circ)
    fast_circ = tofastgates(circ, nqubits)
    return fast_circ
end

"""
    tofastgates(circ::Vector{G}, nqubits::Integer) where {G<:Gate}

Transforms a circuit in the form of a vector of gates to a vector of potentially faster gates where applicable.
"""
function tofastgates(circ::Vector{G}, nqubits::Integer) where {G<:Gate}
    fast_circ = Vector{Gate}(undef, length(circ))
    tofastgates!(fast_circ, nqubits)
    return fast_circ
end

"""
    tofastgates!(circ::Vector{G}) where {G<:Gate}

Transforms a circuit in the form of a vector of gates, converting gates in-place to potentially faster gates where applicable.
The maximum number of qubits is determined from the gates in the circuit, but they all require a `qinds` field.
"""
function tofastgates!(circ::Vector{G}) where {G<:Gate}
    # Find the maximum number of qubits
    nqubits = _getmaxqubits(circ)
    tofastgates!(circ, nqubits)
    return circ
end

"""
    tofastgates!(circ::Vector{G}) where {G<:Gate}

Transforms a circuit in the form of a vector of gates, converting gates in-place to potentially faster gates where applicable.
"""
function tofastgates!(circ::Vector{G}, nqubits::Integer) where {G<:Gate}
    # TODO: This could fail if circ is too concretely typed
    for (ii, gate) in enumerate(circ)
        circ[ii] = tofastgates(gate, nqubits)
    end
    return circ
end


function _getmaxqubits(circ::Vector{G}) where {G<:Gate}
    nqubits = 1
    for gate in circ
        nqubits = max(nqubits, _getmaxqubits(gate))
    end
    return nqubits
end

"""
Get the maximum qubit index that a gate acts on by checking whether it has `qind` or `qinds` fields.
"""
function _getmaxqubits(gate::Gate)
    has_qind_field = hasfield(typeof(gate), :qind)
    has_qinds_field = hasfield(typeof(gate), :qinds)
    if !(has_qind_field || has_qinds_field)
        throw(
            ArgumentError(
                "Gate $(typeof(gate)) does not have `qind` or `qinds` fields defined.
                Use tofastgates!(circ, nqubits) instead."
            )
        )
    end

    if has_qind_field
        return maximum(gate.qind)
    else
        return maximum(gate.qinds)
    end
end