# TODO: We should make these type-stable by using tuples and not vectors
"""
    PauliRotation(symbols::Vector{Symbol}, qinds::Vector{Int})

A parametrized Pauli rotation gate acting on the qubits `qinds` with the Pauli string `symbols`.
"""
struct PauliRotation <: ParametrizedGate
    symbols::Vector{Symbol}
    qinds::Vector{Int}
end

"""
    PauliRotation(symbol::Symbol, qind::Int)

Constructor for a `PauliRotation` acting on a single qubit `qind` with the Pauli `symbol`.
"""
function PauliRotation(symbol::Symbol, qind::Int)
    return PauliRotation([symbol], [qind])
end

"""
    PauliRotation(symbols::Union{AbstractArray,Tuple,Base.Generator}, qinds::Union{AbstractArray,Tuple,Base.Generator})

Constructor for a `PauliRotation` acting on the qubits `qinds` with the Pauli string `symbols`. 
Converts the types of the input arguments to the correct types for `PauliRotation`.
"""
function PauliRotation(symbols::Union{AbstractArray,Tuple,Base.Generator}, qinds::Union{AbstractArray,Tuple,Base.Generator})
    return PauliRotation(collect(symbols), collect(qinds))
end

"""
    PauliRotation(symbols, qinds, theta)

Constructor for a frozen `PauliRotation` acting on the qubits `qinds` with the Pauli string `symbols`, and with fixed parameter `theta`.
"""
function PauliRotation(symbols, qinds, theta)
    return FrozenGate(PauliRotation(symbols, qinds), theta)
end

"""
    MaskedPauliRotation(symbols::Vector{Symbol}, qinds::Vector{Int}, term::PauliStringType)

A parametrized Pauli rotation gate acting on the qubits `qinds` with the Pauli string `symbols`.
The `term` is the integer representation of the Pauli string with the correct integer type for the total number of qubits.
This allows for faster application of the gate.
See `tomaskedpaulirotation` for conversion from `PauliRotation`, which is the easiest way to construct a `MaskedPauliRotation`.
"""
struct MaskedPauliRotation{T} <: ParametrizedGate where {T<:PauliStringType}  # TODO rename
    symbols::Vector{Symbol}
    qinds::Vector{Int}
    generator_mask::T
end

import Base.show
"""
Pretty print for `MaskedPauliRotation`.
"""
function show(io::IO, gate::MaskedPauliRotation)
    print(io, "MaskedPauliRotation{$(typeof(gate.generator_mask))}($(gate.symbols), $(gate.qinds))")
end

"""
    PauliRotationUnion

Union type for `PauliRotation` and `MaskedPauliRotation`, useful for functions which handle either agnostically.
"""
PauliRotationUnion = Union{PauliRotation,MaskedPauliRotation}


"""
    _tomaskedpaulirotation(circ::Vector{Gate})

Returns a circuit where `PauliRotation` gates are transformed to `MaskedPauliRotation` gates.
This allows for significantly faster computation with the gate.
"""
function _tomaskedpaulirotation(circ::Vector{G}, nqubits::Integer) where {G<:Gate}
    TT = getinttype(nqubits)
    return _tomaskedpaulirotation(circ, TT)
end

"""
    _tomaskedpaulirotation(circ::Vector{Gate})

Returns a circuit where `PauliRotation` gates are transformed to `MaskedPauliRotation` gates.
This allows for significantly faster computation with the gate.
"""
function _tomaskedpaulirotation(circ::Vector{G}, ::Type{TT}) where {G<:Gate,TT<:PauliStringType}
    masked_circ = copy(circ)
    for (ii, gate) in enumerate(masked_circ)
        if isa(gate, PauliRotation)
            masked_circ[ii] = _tomaskedpaulirotation(gate, TT)
        end
    end
    return masked_circ
end

"""
    _tomaskedpaulirotation(pauli_gate::PauliRotation, nqubits::Integer)

Transforms a `PauliRotation` to a `MaskedPauliRotation` which carries the integer representation of the gate generator.
This allows for significantly faster computation with the gate.
"""
function _tomaskedpaulirotation(pauli_gate::PauliRotation, nqubits::Integer)
    pstr_term = symboltoint(nqubits, pauli_gate.symbols, pauli_gate.qinds)
    return MaskedPauliRotation(pauli_gate.symbols, pauli_gate.qinds, pstr_term)
end

"""
    _tomaskedpaulirotation(frozen_gate::FrozenGate, nqubits::Integer)

Transforms a `FrozenGate` with a `PauliRotation` to a `MaskedPauliRotation` 
which carries the integer representation of the gate generator.
"""
function _tomaskedpaulirotation(frozen_gate::FrozenGate, nqubits::Integer)
    return FrozenGate(_tomaskedpaulirotation(frozen_gate.gate, nqubits), frozen_gate.theta)
end

"""
    _tomaskedpaulirotation(fast_pauli_gate::MaskedPauliRotation, args...)

Return the `MaskedPauliRotation` gate as is.
"""
function _tomaskedpaulirotation(masked_pauli_gate::MaskedPauliRotation, args...)
    return masked_pauli_gate
end

"""
    _tomaskedpaulirotation(pauli_gate::PauliRotation, ::TermType)

Transforms a `PauliRotation` to a `MaskedPauliRotation` which carries the integer representation of the gate generator.
This allows for significantly faster computation with the gate.
"""
function _tomaskedpaulirotation(pauli_gate::PauliRotation, ::Type{TT}) where {TT<:PauliStringType}
    pstr_term = symboltoint(TT, pauli_gate.symbols, pauli_gate.qinds)
    return MaskedPauliRotation(pauli_gate.symbols, pauli_gate.qinds, pstr_term)
end


### Apply PauliRotation helpers

"""
    splitapply(gate::PauliRotation, pstr::PauliStringType, theta, coefficient=1.0; kwargs...)

Apply a `PauliRotation` with an angle `theta` and a coefficient `coefficient` to an integer Pauli string,
assuming that the gate does not commute with the Pauli string.
Returns two pairs of (pstr, coefficient) as one tuple.
Currently `kwargs` are passed to `applycos` and `applysin` for the Surrogate.
"""
function splitapply(gate::PauliRotationUnion, pstr::PauliStringType, theta, coefficient=1.0; kwargs...)
    coeff1 = _applycos(coefficient, theta; kwargs...)
    new_pstr, sign = _getnewpaulistring(gate, pstr)
    coeff2 = _applysin(coefficient, theta; sign=sign, kwargs...)

    return pstr, coeff1, new_pstr, coeff2
end

"""
    commutes(gate::PauliRotation, pstr::PauliString)

Check if a `PauliRotation` commutes with a `PauliString`.
"""
function commutes(gate::PauliRotation, pstr::PauliString)
    return commutes(gate, pstr.term)
end

"""
    commutes(gate::PauliRotationUnion, pstr)

Check if a `PauliRotation` commutes with an integer Pauli string.
"""
function commutes(gate::PauliRotation, pstr::PauliStringType)
    return sum(!commutes(gate_sym, getpauli(pstr, qind)) for (qind, gate_sym) in zip(gate.qinds, gate.symbols)) % 2 == 0
end

"""
    commutes(gate::MaskedPauliRotation, pstr::PauliStringType)

Check if a `MaskedPauliRotation` commutes with an integer Pauli string.
"""
function commutes(gate::MaskedPauliRotation, pstr::PauliStringType)
    return commutes(gate.generator_mask, pstr)
end

"""
    _applysin(old_coeff::Number, theta; sign=1, kwargs...)

Multiply a numerical coefficient with sin(theta) * sign.
"""
function _applysin(old_coeff::Number, theta; sign=1, kwargs...)
    return old_coeff * sin(theta) * sign
end

"""
    _applycos(old_coeff::Number, theta; sign=1, kwargs...)

Multiply a numerical coefficient with cos(theta) * sign.
"""
function _applycos(old_coeff::Number, theta; sign=1, kwargs...)
    return old_coeff * cos(theta) * sign
end

# TODO: Simplify the following functions
"""
    _applysin(pth::PathProperties, theta; sign=1, kwargs...)

Multiply sin(theta) * sign to the `coeff` field of a `PathProperties` object.
Increments the `nsins` and `freq` fields by 1 if applicable.
"""
function _applysin(pth::T, theta; sign=1, kwargs...) where {T<:PathProperties}
    fields = fieldnames(T)

    @inline function updateval(val, field)
        if field == :coeff
            # apply sin to the `coeff` field
            return _applysin(val, theta; sign=sign, kwargs...)
        elseif field == :nsins
            return val + 1
        elseif field == :freq
            return val + 1
        else
            return val
        end
    end
    return T((updateval(getfield(pth, field), field) for field in fields)...)
end

"""
    _applycos(pth::PathProperties, theta; sign=1, kwargs...)

Multiply cos(theta) * sign to the `coeff` field of a `PathProperties` object.
Increments the `ncos` and `freq` fields by 1 if applicable.
"""
function _applycos(pth::T, theta; sign=1, kwargs...) where {T<:PathProperties}
    fields = fieldnames(T)

    @inline function updateval(val, field)
        if field == :coeff
            # apply cos to the `coeff` field
            return _applycos(val, theta; sign=sign, kwargs...)
        elseif field == :ncos
            return val + 1
        elseif field == :freq
            return val + 1
        else
            return val
        end
    end
    return T((updateval(getfield(pth, field), field) for field in fields)...)
end

"""
    _getnewpaulistring(gate::PauliRotation, pstr::PauliStringType)

Get the new Pauli string after applying a `PauliRotation` to an integer Pauli string,
as well as the corresponding ±1 coefficient.
"""
function _getnewpaulistring(gate::PauliRotation, pstr::PauliStringType)
    new_pstr = copy(pstr)

    total_sign = 1  # this coefficient will be imaginary
    for (qind, gate_sym) in zip(gate.qinds, gate.symbols)
        sign, new_partial_str = pauliprod(gate_sym, getpauli(new_pstr, qind))
        total_sign *= sign
        new_pstr = setpauli(new_pstr, new_partial_str, qind)
    end
    return new_pstr, real(1im * total_sign)
end

"""
    _getnewpaulistring(gate::MaskedPauliRotation, pstr::PauliStringType)

Get the new Pauli string after applying a `MaskedPauliRotation` to an integer Pauli string,
as well as the corresponding ±1 coefficient.
"""
function _getnewpaulistring(gate::MaskedPauliRotation, pstr::PauliStringType)
    # TODO: This allocates memory
    sign, new_pstr = pauliprod(gate.generator_mask, pstr, gate.qinds)
    return new_pstr, real(1im * sign)
end

function _incrementcosandfreq(coeff)
    return coeff
end

function _incrementsinandfreq(coeff)
    return coeff
end

function _incrementcosandfreq(pth::T) where {T<:PathProperties}
    fields = fieldnames(T)
    @inline function updateval(val, field)
        if field == :ncos
            return val + 1
        elseif field == :freq
            return val + 1
        else
            return val
        end
    end
    return T((updateval(getfield(pth, field), field) for field in fields)...)
end

function _incrementsinandfreq(pth::T) where {T<:PathProperties}
    fields = fieldnames(T)
    @inline function updateval(val, field)
        if field == :nsins
            return val + 1
        elseif field == :freq
            return val + 1
        else
            return val
        end
    end
    return T((updateval(getfield(pth, field), field) for field in fields)...)
end