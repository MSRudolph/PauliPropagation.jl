struct PauliGate <: ParametrizedGate
    symbols::Vector{Symbol}
    qinds::Vector{Int}
end

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

function commutes(gate::PauliGateUnion, oper)
    return sum(!commutes(gate_sym, getpaulielement(oper, qind)) for (qind, gate_sym) in zip(gate.qinds, gate.symbols)) % 2 == 0
end

function commutes(gate::FastPauliGate, oper::Integer)
    return commutes(gate.bitoperator, oper)
end

function applysin(old_coeff::Number, theta; sign=1, kwargs...)
    return old_coeff * sin(theta) * sign
end

function applycos(old_coeff::Number, theta; sign=1, kwargs...)
    return old_coeff * cos(theta) * sign
end

function applysin(path_properties::PathProperties, theta; sign=1, kwargs...)
    # path_properties = copy(path_properties) # copy not necesasry. Was done in applycos.
    path_properties.nsins += 1
    path_properties.freq += 1

    path_properties.coeff = applysin(path_properties.coeff, theta; sign, kwargs...)
    return path_properties
end

function applycos(path_properties::PathProperties, theta; sign=1, kwargs...)
    path_properties = copy(path_properties)
    path_properties.ncos += 1
    path_properties.freq += 1

    path_properties.coeff = applycos(path_properties.coeff, theta; sign, kwargs...)
    return path_properties
end


function getnewoperator(gate::PauliGate, oper)
    new_oper = copy(oper)

    total_sign = 1  # this coefficient will be imaginary
    for (qind, gate_sym) in zip(gate.qinds, gate.symbols)
        sign, new_partial_op = pauliprod(gate_sym, getpaulielement(new_oper, qind))
        total_sign *= sign
        new_oper = setpaulielement!(new_oper, qind, new_partial_op)
    end
    return real(1im * total_sign), new_oper
end

function getnewoperator(gate::FastPauliGate, oper)
    sign, new_op = pauliprod(gate.bitoperator, oper, gate.qinds)
    return real(1im * sign), new_op
end