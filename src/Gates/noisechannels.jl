# Depolarzing noise channel

struct DepolarizingNoise <: ParametrizedGate
    qinds::Vector{Int}
end

function apply(gate::DepolarizingNoise, operator, p, coefficient=1.0)

    for qind in gate.qinds
        if getelement(operator, qind) != 0   # non-identity operator
            coefficient *= (1 - p)
        end
    end

    return operator, coefficient
end



### Pauli noise channels

struct PauliXNoise <: ParametrizedGate
    qinds::Vector{Int}
end

function apply(gate::PauliXNoise, operator, p, coefficient=1.0)

    for qind in gate.qinds
        if getelement(operator, qind) == 1   # X operator
            coefficient *= (1 - p)
        end
    end

    return operator, coefficient
end


struct PauliYNoise <: ParametrizedGate
    qinds::Vector{Int}
end

function apply(gate::PauliYNoise, operator, p, coefficient=1.0)

    for qind in gate.qinds
        if getelement(operator, qind) == 2   # Y operator
            coefficient *= (1 - p)
        end
    end

    return operator, coefficient
end


struct PauliZNoise <: ParametrizedGate
    qinds::Vector{Int}
end

function apply(gate::PauliZNoise, operator, p, coefficient=1.0)

    for qind in gate.qinds
        if getelement(operator, qind) == 3   # Z operator
            coefficient *= (1 - p)
        end
    end

    return operator, coefficient
end

