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


"""
    actsdiagonally(gate::AmplitudeDampingNoise, pstr::PauliStringType)

Check if the amplitude damping noise channel acts diagonally on the Pauli string `pstr`.
This implies no splitting (acting diagonally) which happens when acting on I, X, and Y.
"""
function actsdiagonally(gate::AmplitudeDampingNoise, pstr::PauliStringType; kwargs...)
    return getpauli(pstr, gate.qind) != 3
end

"""
    diagonalapply(gate::AmplitudeDampingNoise, pstr::PauliStringType, gamma, coefficient=1.0)

Apply an amplitude damping noise channel to an integer Pauli string `pstr` with noise strength `gamma`.
This is under the assumption that it has been checked that the noise channel acts diagonally on the Pauli string.
Returns a tuple of Pauli string and coefficient.
Physically `gamma` is restricted to the range `[0, 1]`.
A coefficient of the Pauli string can optionally be passed as `coefficient`.
"""
function diagonalapply(gate::AmplitudeDampingNoise, pstr::PauliStringType, gamma, coefficient=1.0; kwargs...)

    local_pauli = getpauli(pstr, gate.qind)

    if local_pauli != 0  # non-identity Pauli
        coefficient *= sqrt(1 - gamma)
    end

    return pstr, coefficient
end

"""
    splitapply(gate::AmplitudeDampingNoise, pstr::PauliStringType, gamma, coefficient=1.0)

Apply an amplitude damping noise channel to an integer Pauli string `pstr` with noise strength `gamma`.
This is under the assumption that it has been checked that the noise channel acts on a Z Pauli and splits.
Returns a tuple of two pairs of Pauli strings and coefficients.
Physically `gamma` is restricted to the range `[0, 1]`.
A coefficient of the Pauli string can optionally be passed as `coefficient`.
"""
function splitapply(gate::AmplitudeDampingNoise, pstr::PauliStringType, gamma, coefficient=1.0; kwargs...)
    new_pstr = setpauli(pstr, 0, gate.qind)
    return pstr, (1 - gamma) * coefficient, new_pstr, gamma * coefficient
end