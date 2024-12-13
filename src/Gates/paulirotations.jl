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
    FastPauliRotation(symbols::Vector{Symbol}, qinds::Vector{Int}, term::PauliStringType)

A parametrized Pauli rotation gate acting on the qubits `qinds` with the Pauli string `symbols`.
The `term` is the integer representation of the Pauli string with the correct integer type for the total number of qubits.
This allows for faster application of the gate.
See `tofastpaulirotation` for conversion from `PauliRotation`, which is the easiest way to construct a `FastPauliRotation`.
"""
struct FastPauliRotation{T} <: ParametrizedGate where {T<:PauliStringType}  # TODO rename
    symbols::Vector{Symbol}
    qinds::Vector{Int}
    term::T
end

import Base.show
"""
Pretty print for `FastPauliRotation`.
"""
function show(io::IO, gate::FastPauliRotation)
    print(io, "FastPauliRotation{$(typeof(gate.term))}($(gate.symbols), $(gate.qinds))")
end

"""
    PauliRotationUnion

Union type for `PauliRotation` and `FastPauliRotation`, useful for functions which handle either agnostically.
"""
PauliRotationUnion = Union{PauliRotation,FastPauliRotation}


"""
    tofastpaulirotation(circ::Vector{Gate})

Returns a circuit where `PauliRotation` gates are transformed to `FastPauliRotation` gates.
This allows for significantly faster computation with the gate.
"""
function tofastpaulirotation(circ::Vector{G}, nqubits::Integer) where {G<:Gate}
    TT = getinttype(nqubits)
    return tofastpaulirotation(circ, TT)
end

"""
    tofastpaulirotation(circ::Vector{Gate})

Returns a circuit where `PauliRotation` gates are transformed to `FastPauliRotation` gates.
This allows for significantly faster computation with the gate.
"""
function tofastpaulirotation(circ::Vector{G}, ::Type{TT}) where {G<:Gate,TT<:PauliStringType}
    fast_circ = copy(circ)
    for (ii, gate) in enumerate(fast_circ)
        if isa(gate, PauliRotation)
            fast_circ[ii] = tofastpaulirotation(gate, TT)
        end
    end
    return fast_circ
end

"""
    tofastpaulirotation(pauli_gate::PauliRotation, nqubits::Integer)

Transforms a `PauliRotation` to a `FastPauliRotation` which carries the integer representation of the gate generator.
This allows for significantly faster computation with the gate.
"""
function tofastpaulirotation(pauli_gate::PauliRotation, nqubits::Integer)
    pstr_term = symboltoint(nqubits, pauli_gate.symbols, pauli_gate.qinds)
    return FastPauliRotation(pauli_gate.symbols, pauli_gate.qinds, pstr_term)
end

"""
    tofastpaulirotation(fast_pauli_gate::PauliRotation, args...)

Return the `FastPauliRotation` gate as is.
"""
function tofastpaulirotation(fast_pauli_gate::FastPauliRotation, args...)
    return fast_pauli_gate
end

"""
    tofastpaulirotation(pauli_gate::PauliRotation, ::TermType)

Transforms a `PauliRotation` to a `FastPauliRotation` which carries the integer representation of the gate generator.
This allows for significantly faster computation with the gate.
"""
function tofastpaulirotation(pauli_gate::PauliRotation, ::Type{TT}) where {TT<:PauliStringType}
    pstr_term = symboltoint(TT, pauli_gate.symbols, pauli_gate.qinds)
    return FastPauliRotation(pauli_gate.symbols, pauli_gate.qinds, pstr_term)
end


### Apply Pauli gates  
# TODO: Implement an apply wrapper function for `PauliString` that works for every low-level apply function.
# like function apply(gate::Gate, pstr::Paulistring, theta, coefficient) and loops over the outcomes to create `PauliString`s
"""
    apply(gate::PauliRotationUnion, pstr::PauliString, theta)

Apply a `(Fast)PauliRotation` with an angle `theta` to a `PauliString`.
Returns either a single `PauliString` or a tuple of two `PauliString`s. 
The latter is the case when the `gate` does not commute with the `PauliString`.
"""
function apply(gate::PauliRotationUnion, pstr::PauliString, theta)
    if commutes(gate, pstr)
        return pstr
    else
        pstr1, c1, pstr2, c2 = applynoncummuting(gate, pstr.term, theta, pstr.coeff)
        return PauliString(pstr.nqubits, pstr1, c1), PauliString(pstr.nqubits, pstr2, c2)
    end
end

"""
    apply(gate::PauliRotationUnion, pstr::PauliStringType, theta, coefficient=1.0)

Apply a `(Fast)PauliRotation` with an angle `theta` and a coefficient `coefficient` to an integer Pauli string.
Returns either one pair of (pstr, coefficient) in one tuple or two pairs as one tuple.
The latter is the case when the `gate` does not commute with the Pauli string.
"""
function apply(gate::PauliRotationUnion, pstr::PauliStringType, theta, coefficient=1.0)
    if commutes(gate, pstr)
        return pstr, coefficient
    else
        return applynoncummuting(gate, pstr, theta, coefficient)
    end
end

"""
    applynoncummuting(gate::PauliRotationUnion, pstr::PauliStringType, theta, coefficient=1.0; kwargs...)

Apply a `(Fast)PauliRotation` with an angle `theta` and a coefficient `coefficient` to an integer Pauli string,
assuming that the gate does not commute with the Pauli string.
Returns two pairs of (pstr, coefficient) as one tuple.
Currently `kwargs` are passed to `applycos` and `applysin` for the Surrogate.
"""
function applynoncummuting(gate::PauliRotationUnion, pstr::PauliStringType, theta, coefficient=1.0; kwargs...)
    coeff1 = applycos(coefficient, theta; kwargs...)
    new_pstr, sign = getnewpaulistring(gate, pstr)
    coeff2 = applysin(coefficient, theta; sign=sign, kwargs...)

    return pstr, coeff1, new_pstr, coeff2
end

"""
    commutes(gate::PauliRotationUnion, pstr::PauliString)

Check if a `(Fast)PauliRotation` commutes with a `PauliString`.
"""
function commutes(gate::PauliRotationUnion, pstr::PauliString)
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
    commutes(gate::FastPauliRotation, pstr::PauliStringType)

Check if a `FastPauliRotation` commutes with an integer Pauli string.
"""
function commutes(gate::FastPauliRotation, pstr::PauliStringType)
    return commutes(gate.term, pstr)
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
    getnewpaulistring(gate::PauliRotation, pstr::PauliStringType)

Get the new Pauli string after applying a `PauliRotation` to an integer Pauli string,
as well as the corresponding ±1 coefficient.
"""
function getnewpaulistring(gate::PauliRotation, pstr::PauliStringType)
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
    getnewpaulistring(gate::FastPauliRotation, pstr::PauliStringType)

Get the new Pauli string after applying a `FastPauliRotation` to an integer Pauli string,
as well as the corresponding ±1 coefficient.
"""
function getnewpaulistring(gate::FastPauliRotation, pstr::PauliStringType)
    # TODO: This allocates memory
    sign, new_pstr = pauliprod(gate.term, pstr, gate.qinds)
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