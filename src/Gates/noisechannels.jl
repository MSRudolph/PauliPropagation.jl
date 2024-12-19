### noisechannels.jl
##
# A file for noise channels. 
# In particular Pauli noise channels and amplitude damping noise.
##
###


# Depolarzing noise channel
"""
Abstract type for parametrized noise channels
"""
abstract type ParametrizedNoiseChannel <: ParametrizedGate end

"""
Abstract type for Pauli noise, i.e., noise that is diagonal in Pauli basis
"""
abstract type PauliNoise <: ParametrizedNoiseChannel end


"""
    DepolarizingNoise(qind::Int)

A depolarizing noise channel acting on the qubit at index `qind`.
Will damp X, Y, and Z Paulis equally.
"""
struct DepolarizingNoise <: PauliNoise
    qind::Int
end

"""
    DepolarizingNoise(qind::Int, p::Real)

A frozen depolarizing noise channel acting on the qubit at index `qind` with noise strength `p`.
Will damp X, Y, and Z Paulis equally.
"""
function DepolarizingNoise(qind::Int, p::Real)
    return FrozenGate(DepolarizingNoise(qind), p)
end

function isdamped(::DepolarizingNoise, pauli::PauliType)
    return pauli != 0
end

"""
    DephasingNoise(qind::Int)

A dephasing noise channel acting on the qubit at index `qind`.
Will damp X and Y Paulis equally.
"""
struct DephasingNoise <: PauliNoise
    qind::Int
end

"""
    DephasingNoise(qind::Int, p::Real)

A frozen dephasing noise channel acting on the qubit at index `qind` with noise strength `p`.
Will damp X and Y Paulis equally.
"""
function DephasingNoise(qind::Int, p::Real)
    return FrozenGate(DephasingNoise(qind), p)
end

function isdamped(::DephasingNoise, pauli::PauliType)
    return pauli == 1 || pauli == 2
end

### Individual Pauli noise channels
"""
    PauliXNoise(qind::Int)

A Pauli-X noise channel acting on the qubit at index `qind`.
Will damp X Paulis.
"""
struct PauliXNoise <: PauliNoise
    qind::Int
end

"""
    PauliXNoise(qind::Int, p::Real)

A frozen Pauli-X noise channel acting on the qubit at index `qind` with noise strength `p`.
Will damp X Paulis.
"""
function PauliXNoise(qind::Int, p::Real)
    return FrozenGate(PauliXNoise(qind), p)
end

function isdamped(::PauliXNoise, pauli::PauliType)
    return pauli == 1
end

"""
    PauliYNoise(qind::Int)

A Pauli-Y noise channel acting on the qubit at index `qind`.
Will damp Y Paulis.
"""
struct PauliYNoise <: PauliNoise
    qind::Int
end

"""
    PauliYNoise(qind::Int, p::Real)

A frozen Pauli-Y noise channel acting on the qubit at index `qind` with noise strength `p`.
Will damp Y Paulis.
"""
function PauliYNoise(qind::Int, p::Real)
    return FrozenGate(PauliYNoise(qind), p)
end

function isdamped(::PauliYNoise, pauli::PauliType)
    return pauli == 2
end

"""
    PauliZNoise(qind::Int)

A Pauli-Z noise channel acting on the qubit at index `qind`.
Will damp Z Paulis.
"""
struct PauliZNoise <: PauliNoise
    qind::Int
end

"""
    PauliZNoise(qind::Int, p::Real)

A frozen Pauli-Z noise channel acting on the qubit at index `qind` with noise strength `p`.
Will damp Z Paulis.
"""
function PauliZNoise(qind::Int, p::Real)
    return FrozenGate(PauliZNoise(qind), p)
end

function isdamped(::PauliZNoise, pauli::PauliType)
    return pauli == 3
end

"""
    AmplitudeDampingNoise(qind::Int)

An amplitude damping noise channel acting on the qubit at index `qind`.
Damps X and Y Paulis, and splits Z into and I and Z component (in the transposed Heisenberg picture).
"""
struct AmplitudeDampingNoise <: ParametrizedGate
    qind::Int
end

"""
    AmplitudeDampingNoise(qind::Int, gamma::Real)

A frozen amplitude damping noise channel acting on the qubit at index `qind` with noise strength `gamma`.
Damps X and Y Paulis, and splits Z into and I and Z component (in the transposed Heisenberg picture).
"""
function AmplitudeDampingNoise(qind::Int, gamma::Real)
    return FrozenGate(AmplitudeDampingNoise(qind), gamma)
end
