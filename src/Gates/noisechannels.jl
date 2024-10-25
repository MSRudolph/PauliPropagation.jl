# Depolarzing noise channel

abstract type PauliNoise <: ParametrizedGate end

struct DepolarizingNoise <: PauliNoise
    qind::Int
end

function apply(gate::DepolarizingNoise, operator, p, coefficient=1.0)

    if getpaulielement(operator, gate.qind) != 0   # non-identity operator
        coefficient *= (1 - p)
    end

    return operator, coefficient
end


### Pauli noise channels

struct PauliXNoise <: PauliNoise
    qind::Int
end

function apply(gate::PauliXNoise, operator, p, coefficient=1.0)

    if getpaulielement(operator, gate.qind) == 1   # X operator
        coefficient *= (1 - p)
    end

    return operator, coefficient
end


struct PauliYNoise <: PauliNoise
    qind::Int
end

function apply(gate::PauliYNoise, operator, p, coefficient=1.0)

    if getpaulielement(operator, gate.qind) == 2   # Y operator
        coefficient *= (1 - p)
    end

    return operator, coefficient
end


struct PauliZNoise <: PauliNoise
    qind::Int
end

function apply(gate::PauliZNoise, operator, p, coefficient=1.0)

    if getpaulielement(operator, gate.qind) == 3   # Z operator
        coefficient *= (1 - p)
    end

    return operator, coefficient
end


struct AmplitudeDampingNoise <: ParametrizedGate
    qind::Int
end

function apply(gate::AmplitudeDampingNoise, operator, p, coefficient=1.0)

    if actsdiagonally(gate, operator)  # test for Z operator
        return diagonalapply(gate, operator, p, coefficient)
    end

    return splitapply(gate, operator, p, coefficient)
end

function actsdiagonally(gate::AmplitudeDampingNoise, operator)
    return getpaulielement(operator, gate.qind) != 3
end

function diagonalapply(gate::AmplitudeDampingNoise, operator, p, coefficient=1.0)

    local_op = getpaulielement(operator, gate.qind)

    if local_op != 0  # non-identity operator
        coefficient *= sqrt(1 - p)
    end

    return operator, coefficient
end

function splitapply(gate::AmplitudeDampingNoise, operator, p, coefficient=1.0)
    new_operator = copy(operator)
    new_operator = setpaulielement!(new_operator, gate.qind, 0)
    return operator, (1 - p) * coefficient, new_operator, p * coefficient
end