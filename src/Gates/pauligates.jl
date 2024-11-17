# TODO: We should make these type-stable by using tuples and not vectors
"""
    PauliGate(symbols::Vector{Symbol}, qinds::Vector{Int})

A parametrized Pauli rotation gate acting on the qubits `qinds` with the Pauli operators `symbols`.
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

Constructor for a `PauliGate` acting on the qubits `qinds` with the Pauli operators `symbols`. 
Converts the types of the input arguments to the correct types for `PauliGate`.
"""
function PauliGate(symbols::Union{AbstractArray,Tuple,Base.Generator}, qinds::Union{AbstractArray,Tuple,Base.Generator})
    # `collect` returns a vector for any iterable 
    return PauliGate(collect(symbols), collect(qinds))
end

"""
    FastPauliGate(symbols::Vector{Symbol}, qinds::Vector{Int}, bitoperator::PauliStringType)

A parametrized Pauli rotation gate acting on the qubits `qinds` with the Pauli operators `symbols`.
The `bitoperator` is the integer representation of the Pauli operators with the correct integer type for the total number of qubits.
This allows for faster application of the gate.
See `tofastgates` for conversion from `PauliGate`, which is the easiest way to construct a `FastPauliGate`.
"""
struct FastPauliGate{T} <: ParametrizedGate where {T<:PauliStringType}
    symbols::Vector{Symbol}
    qinds::Vector{Int}
    bitoperator::T
end

import Base.show
function show(io::IO, fastpauligate::FastPauliGate)
    print(io, "FastPauliGate{$(typeof(fastpauligate.bitoperator))}($(fastpauligate.symbols), $(fastpauligate.qinds))")
end

"""
    PauliGateUnion

Union type for `PauliGate` and `FastPauliGate`.
"""
PauliGateUnion = Union{PauliGate,FastPauliGate}

"""
    tofastgates(gate::Gate, nqubits::Integer)

Transforms a gate to a potentially faster but more involved gate type. 
This is not implemented for all gates.
"""
tofastgates(gate::Gate, nqubits::Integer) = gate  # TODO: move this to a more general place

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

"""
    tofastgates(circ::Vector{G}) where {G<:Gate}

Transforms a circuit in the form of a vector of gates to a vector of potentially faster gates where applicable.
"""
function tofastgates(circ::Vector{G}) where {G<:Gate}
    # Find the maximum number of qubits
    nq = 1
    for gate in circ
        nq = max(nq, maximum(gate.qinds))
    end

    fast_circ = similar(circ)
    for (ii, gate) in enumerate(circ)
        fast_circ[ii] = tofastgates(gate, nq)
    end
    return fast_circ
end

"""
    tofastgates!(circ::Vector{G}) where {G<:Gate}

Transforms a circuit in the form of a vector of gates, converting gates in-place to potentially faster gates where applicable.
"""
function tofastgates!(circ::Vector{G}) where {G<:Gate}
    # Find the maximum number of qubits
    nq = 1
    for gate in circ
        nq = max(nq, maximum(gate.qinds))
    end

    # TODO: This could fail if circ is too concretely typed
    for (ii, gate) in enumerate(circ)
        circ[ii] = tofastgates(gate, nq)
    end
    return circ
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
    new_pstr, sign = getnewoperator(gate, pstr)
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

"""
    applysin(path_properties::PathProperties, theta; sign=1, kwargs...)

Multiply sin(theta) * sign to the coefficient of a `PathProperties` object.
Increments the `nsins` and `freq` fields by 1.
"""
function applysin(path_properties::PathProperties, theta; sign=1, kwargs...)
    # path_properties = copy(path_properties) # copy not necessary. Was done in applycos.
    _incrementsinandfreq!(path_properties)

    path_properties.coeff = applysin(path_properties.coeff, theta; sign, kwargs...)
    return path_properties
end

"""
    applycos(path_properties::PathProperties, theta; sign=1, kwargs...)

Multiply cos(theta) * sign to the coefficient of a `PathProperties` object.
Increments the `ncos` and `freq` fields by 1.
"""
function applycos(path_properties::PathProperties, theta; sign=1, kwargs...)
    path_properties = copy(path_properties)
    _incrementcosandfreq!(path_properties)

    path_properties.coeff = applycos(path_properties.coeff, theta; sign, kwargs...)
    return path_properties
end

"""
    getnewoperator(gate::PauliGate, pstr::PauliStringType)

Get the new Pauli string after applying a `PauliGate` to an integer Pauli string,
as well as the corresponding ±1 coefficient.
"""
function getnewoperator(gate::PauliGate, pstr::PauliStringType)
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
    getnewoperator(gate::FastPauliGate, pstr::PauliStringType)

Get the new Pauli string after applying a `FastPauliGate` to an integer Pauli string,
as well as the corresponding ±1 coefficient.
"""
function getnewoperator(gate::FastPauliGate, pstr::PauliStringType)
    # TODO: This allocates memory
    sign, new_pstr = pauliprod(gate.bitoperator, pstr, gate.qinds)
    return new_pstr, real(1im * sign)
end

function _incrementcosandfreq!(coeff)
    return
end

function _incrementsinandfreq!(coeff)
    return
end

function _incrementcosandfreq!(path_properties::PathProperties)
    path_properties.ncos += 1
    path_properties.freq += 1
    return
end

function _incrementsinandfreq!(path_properties::PathProperties)
    path_properties.nsins += 1
    path_properties.freq += 1
    return
end