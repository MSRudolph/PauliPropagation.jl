### fastgates.jl
##
# This is an interface to convert gates to potentially more performant versions of that gate.
# Often the faster gates are not as nice to handle at a high level, so we allow for conversion.
##
###

"""
    tofastgates(gate::Gate, nqubits::Integer)

Transforms a gate to a potentially faster but more involved gate type. 
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
    fast_circ = deepcopy(circ)
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

function _maxqubits(gate::G) where {G<:Gate}
    # check for typical qind/qinds fields
    if hasfield(G, :qinds)
        return maximum(gate.qinds)
    elseif hasfield(G, :qind)
        return gate.qind
    else
        return -1  # to show that it failed
    end
end

function _getmaxqubits(circ::Vector{G}) where {G<:Gate}
    nqubits = 1
    for gate in circ
        maxqubits = _maxqubits(gate)
        if maxqubits == -1
            throw(
                ArgumentError(
                    "Gate $(typeof(gate)) does not have `qind` or `qinds` field defined.
                    Use tofastgates!(circ, nqubits) instead."
                )
            )
        end
        nqubits = max(nqubits, maxqubits)
    end
    return nqubits
end