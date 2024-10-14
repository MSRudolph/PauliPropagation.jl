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

tofastgates(gate::Gate, nq::Integer) = gate

function tofastgates(pauli_gate::PauliGate, nq::Integer)
    base_ops = [:I for _ in 1:nq]
    for (qind, op) in zip(pauli_gate.qinds, pauli_gate.symbols)
        base_ops[qind] = op
    end
    return FastPauliGate(pauli_gate.symbols, pauli_gate.qinds, symboltoint(base_ops))
end

function tofastgates(circ::AbstractVector)
    nq = 1
    for gate in circ
        nq = max(nq, maximum(gate.qinds))
    end
    fast_circ = Gate[]

    for gate in circ
        push!(fast_circ, tofastgates(gate, nq))
    end
    return fast_circ
end


### Apply Pauli gates  

function apply(gate::PauliGateUnion, operator, theta, coefficient=1.0)  # TODO: make order of arguments consistent with other apply functions
    if commutes(gate, operator)
        return operator, coefficient
    else
        return applynoncummuting(gate, operator, theta, coefficient)
    end
end

function applynoncummuting(gate::PauliGateUnion, operator, theta, coefficient=1.0; param_idx=0)
    coeff1 = applycos(coefficient, theta; param_idx=param_idx)
    sign, new_oper = getnewoperator(gate, operator)
    coeff2 = applysin(coefficient, theta; sign=sign, param_idx=param_idx)

    return operator, coeff1, new_oper, coeff2   # TODO: when does this ever not allocate memory?
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