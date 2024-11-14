struct PauliGate <: ParametrizedGate
    symbols::Vector{Symbol}
    qinds::Vector{Int}
end

# TODO: We should make these type-stable by using tuples and not vectors
struct FastPauliGate{T} <: ParametrizedGate where {T<:Integer}
    symbols::Vector{Symbol}
    qinds::Vector{Int}
    bitoperator::T
end

PauliGateUnion = Union{PauliGate,FastPauliGate}

tofastgates(gate::Gate, nq::Integer) = gate  # TODO: move this to a more general place

function tofastgates(pauli_gate::PauliGate, nq::Integer)
    base_ops = [:I for _ in 1:nq]
    for (qind, op) in zip(pauli_gate.qinds, pauli_gate.symbols)
        base_ops[qind] = op
    end
    return FastPauliGate(pauli_gate.symbols, pauli_gate.qinds, symboltoint(base_ops))
end

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
function apply(gate::PauliGateUnion, pstr::PauliString, theta)
    if commutes(gate, pstr)
        return pstr
    else
        pstr1, c1, pstr2, c2 = applynoncummuting(gate, pstr.operator, theta, pstr.coeff)
        return PauliString(pstr.nqubits, pstr1, c1), PauliString(pstr.nqubits, pstr2, c2)
    end
end

function apply(gate::PauliGateUnion, operator, theta, coefficient=1.0)
    if commutes(gate, operator)
        return operator, coefficient
    else
        return applynoncummuting(gate, operator, theta, coefficient)
    end
end

function applynoncummuting(gate::PauliGateUnion, operator, theta, coefficient=1.0; kwargs...)
    coeff1 = applycos(coefficient, theta; kwargs...)
    sign, new_oper = getnewoperator(gate, operator)
    coeff2 = applysin(coefficient, theta; sign=sign, kwargs...)

    return operator, coeff1, new_oper, coeff2
end

function commutes(gate::PauliGateUnion, pstr::PauliString)
    return commutes(gate, pstr.operator)
end

function commutes(gate::PauliGateUnion, pstr)
    return sum(!commutes(gate_sym, getpauli(pstr, qind)) for (qind, gate_sym) in zip(gate.qinds, gate.symbols)) % 2 == 0
end

function commutes(gate::FastPauliGate, pstr::PauliStringType)
    return commutes(gate.bitoperator, pstr)
end

function applysin(old_coeff::Number, theta; sign=1, kwargs...)
    return old_coeff * sin(theta) * sign
end

function applycos(old_coeff::Number, theta; sign=1, kwargs...)
    return old_coeff * cos(theta) * sign
end

function applysin(path_properties::PathProperties, theta; sign=1, kwargs...)
    # path_properties = copy(path_properties) # copy not necesasry. Was done in applycos.
    _incrementsinandfreq!(path_properties)

    path_properties.coeff = applysin(path_properties.coeff, theta; sign, kwargs...)
    return path_properties
end

function applycos(path_properties::PathProperties, theta; sign=1, kwargs...)
    path_properties = copy(path_properties)
    _incrementcosandfreq!(path_properties)

    path_properties.coeff = applycos(path_properties.coeff, theta; sign, kwargs...)
    return path_properties
end


function getnewoperator(gate::PauliGate, pstr)
    new_oper = copy(pstr)

    total_sign = 1  # this coefficient will be imaginary
    for (qind, gate_sym) in zip(gate.qinds, gate.symbols)
        sign, new_partial_op = pauliprod(gate_sym, getpauli(new_oper, qind))
        total_sign *= sign
        new_oper = setpauli(new_oper, new_partial_op, qind)
    end
    return real(1im * total_sign), new_oper
end

function getnewoperator(gate::FastPauliGate, oper)
    # TODO: This allocates memory
    sign, new_op = pauliprod(gate.bitoperator, oper, gate.qinds)
    return real(1im * sign), new_op
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