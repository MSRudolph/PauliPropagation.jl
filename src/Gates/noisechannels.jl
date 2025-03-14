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



struct DepolarizingNoise <: PauliNoise
    qind::Int

    @doc """
        DepolarizingNoise(qind::Int)

    A depolarizing noise channel acting on the qubit at index `qind`.
    Will damp X, Y, and Z Paulis equally by a factor of `1-p`.
    """
    DepolarizingNoise(qind::Int) = (_qinds_check(qind); new(qind))
end

"""
    DepolarizingNoise(qind::Int, p::Real)

A frozen depolarizing noise channel acting on the qubit at index `qind` with noise strength `p`.
Will damp X, Y, and Z Paulis equally by a factor of `1-p`.
"""
function DepolarizingNoise(qind::Int, p::Real)
    _check_noise_strength(DepolarizingNoise, p)

    return FrozenGate(DepolarizingNoise(qind), p)
end

function isdamped(::DepolarizingNoise, pauli::PauliType)
    return pauli != 0
end


struct DephasingNoise <: PauliNoise
    qind::Int

    @doc """
        DephasingNoise(qind::Int)

    A dephasing noise channel acting on the qubit at index `qind`.
    Will damp X and Y Paulis equally by a factor of `1-p`.
    """
    DephasingNoise(qind::Int) = (_qinds_check(qind); new(qind))
end

"""
    DephasingNoise(qind::Int, p::Real)

A frozen dephasing noise channel acting on the qubit at index `qind` with noise strength `p`.
Will damp X and Y Paulis equally by a factor of `1-p`.
"""
function DephasingNoise(qind::Int, p::Real)
    _check_noise_strength(DephasingNoise, p)

    return FrozenGate(DephasingNoise(qind), p)
end

function isdamped(::DephasingNoise, pauli::PauliType)
    return pauli == 1 || pauli == 2
end

### Individual Pauli noise damping

struct PauliXDamping <: PauliNoise
    qind::Int

    @doc """
        PauliXDamping(qind::Int)

    A Pauli-X noise damping acting on the qubit at index `qind`.
    Will damp X Paulis by a factor of `1-p`. This alone is not a valid quantum channel.
    """
    PauliXDamping(qind::Int) = (_qinds_check(qind); new(qind))
end

"""
    PauliXDamping(qind::Int, p::Real)

A frozen Pauli-X noise damping acting on the qubit at index `qind` with noise strength `p`.
Will damp X Paulis by a factor of `1-p`. This alone is not a valid quantum channel.
"""
function PauliXDamping(qind::Int, p::Real)
    _check_noise_strength(PauliXDamping, p)

    return FrozenGate(PauliXDamping(qind), p)
end

function isdamped(::PauliXDamping, pauli::PauliType)
    return pauli == 1
end


struct PauliYDamping <: PauliNoise
    qind::Int

    @doc """
        PauliYDamping(qind::Int)

    A Pauli-Y noise damping acting on the qubit at index `qind`.
    Will damp Y Paulis by a factor of `1-p`. This alone is not a valid quantum channel.
    """
    PauliYDamping(qind::Int) = (_qinds_check(qind); new(qind))
end

"""
    PauliYDamping(qind::Int, p::Real)

A frozen Pauli-Y damping acting on the qubit at index `qind` with noise strength `p`.
Will damp Y Paulis by a factor of `1-p`. This alone is not a valid quantum channel.
"""
function PauliYDamping(qind::Int, p::Real)
    _check_noise_strength(PauliYDamping, p)

    return FrozenGate(PauliYDamping(qind), p)
end

function isdamped(::PauliYDamping, pauli::PauliType)
    return pauli == 2
end


struct PauliZDamping <: PauliNoise
    qind::Int

    @doc """
        PauliZDamping(qind::Int)

    A Pauli-Z noise damping acting on the qubit at index `qind`.
    Will damp Z Paulis by a factor of `1-p`. This alone is not a valid quantum channel.
    """
    PauliZDamping(qind::Int) = (_qinds_check(qind); new(qind))
end

"""
    PauliZDamping(qind::Int, p::Real)

A frozen Pauli-Z noise damping acting on the qubit at index `qind` with noise strength `p`.
Will damp Z Paulis by a factor of `1-p`. This alone is not a valid quantum channel.
"""
function PauliZDamping(qind::Int, p::Real)
    _check_noise_strength(PauliZDamping, p)

    return FrozenGate(PauliZDamping(qind), p)
end

function isdamped(::PauliZDamping, pauli::PauliType)
    return pauli == 3
end

"""
    AmplitudeDampingNoise(qind::Int)

An amplitude damping noise channel acting on the qubit at index `qind`.
Damps X and Y Paulis by a factor of sqrt(1-gamma)
and splits Z into and gamma * I and (1-gamma) * Z component (in the transposed Heisenberg picture).
"""
struct AmplitudeDampingNoise <: ParametrizedNoiseChannel
    qind::Int
end

"""
    AmplitudeDampingNoise(qind::Int, gamma::Real)

A frozen amplitude damping noise channel acting on the qubit at index `qind` with noise strength `gamma`.
Damps X and Y Paulis, and splits Z into and I and Z component (in the transposed Heisenberg picture).
"""
function AmplitudeDampingNoise(qind::Int, gamma::Real)
    _check_noise_strength(AmplitudeDampingNoise, gamma)

    return FrozenGate(AmplitudeDampingNoise(qind), gamma)
end


function _check_noise_strength(::Type{G}, p) where {G<:ParametrizedNoiseChannel}
    if !(0 <= p <= 1)
        throw(ArgumentError("$G parameter must be between 0 and 1. Got $p."))
    end
end