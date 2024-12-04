# TODO: We should make these type-stable by using tuples and not vectors
"""
    PauliGate(symbols::Vector{Symbol}, qinds::Vector{Int})

A parametrized Pauli rotation gate acting on the qubits `qinds` with the Pauli string `symbols`.
"""
struct PauliGate <: ParametrizedGate
    symbols::Vector{Symbol}
    qinds::Vector{Int}
end

"""
    PauliGate(symbol::Symbol, qind::Int)

Constructor for a `PauliGate` acting on a single qubit `qind` with the Pauli `symbol`.
"""
function PauliGate(symbol::Symbol, qind::Int)
    return PauliGate([symbol], [qind])
end

"""
    PauliGate(symbols::Union{AbstractArray,Tuple,Base.Generator}, qinds::Union{AbstractArray,Tuple,Base.Generator})

Constructor for a `PauliGate` acting on the qubits `qinds` with the Pauli string `symbols`. 
Converts the types of the input arguments to the correct types for `PauliGate`.
"""
function PauliGate(symbols::Union{AbstractArray,Tuple,Base.Generator}, qinds::Union{AbstractArray,Tuple,Base.Generator})
    return PauliGate(collect(symbols), collect(qinds))
end

"""
    PauliGate(symbols, qinds, theta)

Constructor for a frozen `PauliGate` acting on the qubits `qinds` with the Pauli string `symbols`, and with fixed parameter `theta`.
"""
function PauliGate(symbols, qinds, theta)
    return FrozenGate(PauliGate(symbols, qinds), theta)
end

"""
    FastPauliGate(symbols::Vector{Symbol}, qinds::Vector{Int}, bitoperator::PauliStringType)

A parametrized Pauli rotation gate acting on the qubits `qinds` with the Pauli string `symbols`.
The `bitoperator` is the integer representation of the Pauli string with the correct integer type for the total number of qubits.
This allows for faster application of the gate.
See `tofastgates` for conversion from `PauliGate`, which is the easiest way to construct a `FastPauliGate`.
"""
struct FastPauliGate{T} <: ParametrizedGate where {T<:PauliStringType}  # TODO rename
    symbols::Vector{Symbol}
    qinds::Vector{Int}
    bitoperator::T
end

import Base.show
"""
Pretty print for `FastPauliGate`.
"""
function show(io::IO, fastpauligate::FastPauliGate)
    print(io, "FastPauliGate{$(typeof(fastpauligate.bitoperator))}($(fastpauligate.symbols), $(fastpauligate.qinds))")
end

"""
    PauliGateUnion

Union type for `PauliGate` and `FastPauliGate`.
"""
PauliGateUnion = Union{PauliGate,FastPauliGate}

"""
    tofastgates(pauli_gate::PauliGate, nqubits::Integer)

Transforms a `PauliGate` to a `FastPauliGate` which carries the integer representation of the gate generator.
This allows for significantly faster computation with the gate.
"""
function tofastgates(pauli_gate::PauliGate, nqubits::Integer)
    base_pstr = getinttype(nqubits)(0)
    for (qind, pauli) in zip(pauli_gate.qinds, pauli_gate.symbols)
        base_pstr = setpauli(base_pstr, pauli, qind)
    end
    return FastPauliGate(pauli_gate.symbols, pauli_gate.qinds, base_pstr)
end


### Apply Pauli gates  
# TODO: Implement an apply wrapper function for `PauliString` that works for every low-level apply function.
# like function apply(gate::Gate, pstr::Paulistring, theta, coefficient) and loops over the outcomes to create `PauliString`s
"""
    apply(gate::PauliGateUnion, pstr::PauliString, theta)

Apply a `(Fast)PauliGate` with an angle `theta` to a `PauliString`.
Returns either a single `PauliString` or a tuple of two `PauliString`s. 
The latter is the case when the `gate` does not commute with the `PauliString`.
"""
function apply(gate::PauliGateUnion, pstr::PauliString, theta)
    if commutes(gate, pstr)
        return pstr
    else
        pstr1, c1, pstr2, c2 = applynoncummuting(gate, pstr.operator, theta, pstr.coeff)
        return PauliString(pstr.nqubits, pstr1, c1), PauliString(pstr.nqubits, pstr2, c2)
    end
end

"""
    apply(gate::PauliGateUnion, pstr::PauliStringType, theta, coefficient=1.0)

Apply a `(Fast)PauliGate` with an angle `theta` and a coefficient `coefficient` to an integer Pauli string.
Returns either one pair of (pstr, coefficient) in one tuple or two pairs as one tuple.
The latter is the case when the `gate` does not commute with the Pauli string.
"""
function apply(gate::PauliGateUnion, pstr::PauliStringType, theta, coefficient=1.0)
    if commutes(gate, pstr)
        return pstr, coefficient
    else
        return applynoncummuting(gate, pstr, theta, coefficient)
    end
end

"""
    applynoncummuting(gate::PauliGateUnion, pstr::PauliStringType, theta, coefficient=1.0; kwargs...)

Apply a `(Fast)PauliGate` with an angle `theta` and a coefficient `coefficient` to an integer Pauli string,
assuming that the gate does not commute with the Pauli string.
Returns two pairs of (pstr, coefficient) as one tuple.
Currently `kwargs` are passed to `applycos` and `applysin` for the Surrogate.
"""
function applynoncummuting(gate::PauliGateUnion, pstr::PauliStringType, theta, coefficient=1.0; kwargs...)
    coeff1 = applycos(coefficient, theta; kwargs...)
    new_pstr, sign = getnewpaulistring(gate, pstr)
    coeff2 = applysin(coefficient, theta; sign=sign, kwargs...)

    return pstr, coeff1, new_pstr, coeff2
end

"""
    commutes(gate::PauliGateUnion, pstr::PauliString)

Check if a `(Fast)PauliGate` commutes with a `PauliString`.
"""
function commutes(gate::PauliGateUnion, pstr::PauliString)
    return commutes(gate, pstr.operator)
end

"""
    commutes(gate::PauliGateUnion, pstr)

Check if a `PauliGate` commutes with an integer Pauli string.
"""
function commutes(gate::PauliGate, pstr::PauliStringType)
    return sum(!commutes(gate_sym, getpauli(pstr, qind)) for (qind, gate_sym) in zip(gate.qinds, gate.symbols)) % 2 == 0
end

"""
    commutes(gate::FastPauliGate, pstr::PauliStringType)

Check if a `FastPauliGate` commutes with an integer Pauli string.
"""
function commutes(gate::FastPauliGate, pstr::PauliStringType)
    return commutes(gate.bitoperator, pstr)
end

"""
    applysin(old_coeff::Number, theta; sign=1, kwargs...)

Multiply a numerical coefficient with sin(theta) * sign.
"""
function applysin(old_coeff::Number, theta; sign=1, kwargs...)
    return old_coeff * sin(theta) * sign
end

"""
    applycos(old_coeff::Number, theta; sign=1, kwargs...)

Multiply a numerical coefficient with cos(theta) * sign.
"""
function applycos(old_coeff::Number, theta; sign=1, kwargs...)
    return old_coeff * cos(theta) * sign
end

# TODO: Simplify the following functions
"""
    applysin(pth::PathProperties, theta; sign=1, kwargs...)

Multiply sin(theta) * sign to the `coeff` field of a `PathProperties` object.
Increments the `nsins` and `freq` fields by 1 if applicable.
"""
function applysin(pth::T, theta; sign=1, kwargs...) where {T<:PathProperties}
    fields = fieldnames(T)

    @inline function updateval(val, field)
        if field == :coeff
            # apply sin to the `coeff` field
            return applysin(val, theta; sign=sign, kwargs...)
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
    applycos(pth::PathProperties, theta; sign=1, kwargs...)

Multiply cos(theta) * sign to the `coeff` field of a `PathProperties` object.
Increments the `ncos` and `freq` fields by 1 if applicable.
"""
function applycos(pth::T, theta; sign=1, kwargs...) where {T<:PathProperties}
    fields = fieldnames(T)

    @inline function updateval(val, field)
        if field == :coeff
            # apply cos to the `coeff` field
            return applycos(val, theta; sign=sign, kwargs...)
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
    getnewpaulistring(gate::PauliGate, pstr::PauliStringType)

Get the new Pauli string after applying a `PauliGate` to an integer Pauli string,
as well as the corresponding ±1 coefficient.
"""
function getnewpaulistring(gate::PauliGate, pstr::PauliStringType)
    new_pstr = copy(pstr)

    total_sign = 1  # this coefficient will be imaginary
    for (qind, gate_sym) in zip(gate.qinds, gate.symbols)
        sign, new_partial_op = pauliprod(gate_sym, getpauli(new_pstr, qind))
        total_sign *= sign
        new_pstr = setpauli(new_pstr, new_partial_op, qind)
    end
    return new_pstr, real(1im * total_sign)
end

"""
    getnewpaulistring(gate::FastPauliGate, pstr::PauliStringType)

Get the new Pauli string after applying a `FastPauliGate` to an integer Pauli string,
as well as the corresponding ±1 coefficient.
"""
function getnewpaulistring(gate::FastPauliGate, pstr::PauliStringType)
    # TODO: This allocates memory
    sign, new_pstr = pauliprod(gate.bitoperator, pstr, gate.qinds)
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