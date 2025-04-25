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



## The other Pauli noise channels
# Pauli-X noise channel
struct PauliXNoise <: PauliNoise
    qind::Int

    @doc """
        PauliXNoise(qind::Int)
        PauliXNoise(qind::Int, p::Real)

    A Pauli-X noise channel acting on the qubit at index `qind`.
    If `p` is provided, this returns a frozen gate with that noise strength.
    Will damp Y and Z Paulis equally by a factor of `1-p`.
    This corresponds to inserting a random Pauli X operator into the circuit with probability `p/2`.
    """
    PauliXNoise(qind::Int) = (_qinds_check(qind); new(qind))
end


# the frozen gate version
function PauliXNoise(qind::Int, p::Real)
    _check_noise_strength(DephasingNoise, p)

    return FrozenGate(PauliXNoise(qind), p)
end


function isdamped(::PauliXNoise, pauli::PauliType)
    return pauli == 2 || pauli == 3
end


# Pauli-Y noise channel
struct PauliYNoise <: PauliNoise
    qind::Int

    @doc """
        PauliYNoise(qind::Int)
        PauliYNoise(qind::Int, p::Real)

    A Pauli-Y noise channel acting on the qubit at index `qind`.
    If `p` is provided, this returns a frozen gate with that noise strength.
    Will damp X and Z Paulis equally by a factor of `1-p`.
    This corresponds to inserting a random Pauli Y operator into the circuit with probability `p/2`.
    """
    PauliYNoise(qind::Int) = (_qinds_check(qind); new(qind))
end


# the frozen gate version
function PauliYNoise(qind::Int, p::Real)
    _check_noise_strength(DephasingNoise, p)

    return FrozenGate(PauliYNoise(qind), p)
end


function isdamped(::PauliYNoise, pauli::PauliType)
    return pauli == 1 || pauli == 3
end


# Pauli-Z noise channel
struct PauliZNoise <: PauliNoise
    qind::Int

    @doc """
        PauliZNoise(qind::Int)
        PauliZNoise(qind::Int, p::Real)

    A Pauli-Z noise channel acting on the qubit at index `qind`.
    If `p` is provided, this returns a frozen gate with that noise strength.
    Will damp X and Y Paulis equally by a factor of `1-p`.
    This corresponds to inserting a random Pauli Z operator with probability `p/2`.
    """
    PauliZNoise(qind::Int) = (_qinds_check(qind); new(qind))
end


# the frozen gate version 
function PauliZNoise(qind::Int, p::Real)
    _check_noise_strength(DephasingNoise, p)

    return FrozenGate(PauliZNoise(qind), p)
end


function isdamped(::PauliZNoise, pauli::PauliType)
    return pauli == 1 || pauli == 2
end


## DepolarizingNoise is an alias for PauliZNoise
"""
    DephasingNoise(qind::Int)
    DephasingNoise(qind::Int, p::Real)

This is an alias for `PauliZNoise`.
If `p` is provided, this returns a frozen gate with that noise strength.
A dephasing noise channel acting on the qubit at index `qind`.
Will damp X and Y Paulis equally by a factor of `1-p`.
"""
const DephasingNoise = PauliZNoise


### Individual Pauli noise damping
# these are not exported because they are not valid quantum channels

struct PauliXDamping <: PauliNoise
    qind::Int

    @doc """
        PauliXDamping(qind::Int)
        PauliXDamping(qind::Int, p::Real)

    A Pauli-X noise damping acting on the qubit at index `qind`.
    If `p` is provided, this returns a frozen gate with that damping strength.
    Will damp X Paulis by a factor of `1-p`. 
    This alone is not a valid quantum channel.
    """
    PauliXDamping(qind::Int) = (_qinds_check(qind); new(qind))
end


# the frozen gate version
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
    If `p` is provided, this returns a frozen gate with that damping strength.
    Will damp Y Paulis by a factor of `1-p`. 
    This alone is not a valid quantum channel.
    """
    PauliYDamping(qind::Int) = (_qinds_check(qind); new(qind))
end


# the frozen gate version
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
    If `p` is provided, this returns a frozen gate with that damping strength.
    Will damp Z Paulis by a factor of `1-p`. 
    This alone is not a valid quantum channel.
    """
    PauliZDamping(qind::Int) = (_qinds_check(qind); new(qind))
end


# the frozen gate version
function PauliZDamping(qind::Int, p::Real)
    _check_noise_strength(PauliZDamping, p)

    return FrozenGate(PauliZDamping(qind), p)
end

function isdamped(::PauliZDamping, pauli::PauliType)
    return pauli == 3
end

## Amplitude damping noise
"""
    AmplitudeDampingNoise(qind::Int)
    AmplitudeDampingNoise(qind::Int, gamma::Real)

An amplitude damping noise channel acting on the qubit at index `qind`.
If `gamma` is provided, this returns a frozen gate with that noise strength.
Damps X and Y Paulis by a factor of sqrt(1-gamma)
and splits Z into and gamma * I and (1-gamma) * Z component (in the transposed Heisenberg picture).
"""
struct AmplitudeDampingNoise <: ParametrizedNoiseChannel
    qind::Int
end


# the frozen gate version
function AmplitudeDampingNoise(qind::Int, gamma::Real)
    _check_noise_strength(AmplitudeDampingNoise, gamma)

    return FrozenGate(AmplitudeDampingNoise(qind), gamma)
end


function _check_noise_strength(::Type{G}, p) where {G<:ParametrizedNoiseChannel}
    if !(0 <= p <= 1)
        throw(ArgumentError("$G parameter must be between 0 and 1. Got $p."))
    end
end