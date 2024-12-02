### PAULI GATES
"""
    applygatetoall!(gate::PauliGateUnion, thetas, operator_dict, new_operator_dict, args...; kwargs...)

Overload of `applygatetoall!` for `PauliGate` and `FastPauliGate` gates.
Both `operator_dict` and `new_operator_dict` contain operators which will be merged later.
"""
function applygatetoall!(gate::PauliGateUnion, theta, operator_dict, new_operator_dict, args...; kwargs...)
    # TODO: there is a lot of code duplication. Can we write a more general function?

    for (operator, coeff) in operator_dict
        applygatetoone!(gate, operator, coeff, theta, operator_dict, new_operator_dict; kwargs...)
    end

    return operator_dict, new_operator_dict  # don't swap psums around
end

"""
    applygatetoone!(gate::PauliGateUnion, operator, coefficient, theta, operator_dict, new_operator_dict, args...; kwargs...)

Overload of `applygatetoone!` for `PauliGate` and `FastPauliGate` gates. 
Checks for commutation of `gate` and `operator`, and applies the gate to the operator if they don't.
"""
@inline function applygatetoone!(gate::PauliGateUnion, operator, coefficient, theta, operator_dict, new_operator_dict, args...; kwargs...)

    if commutes(gate, operator)
        return
    end

    operator, coeff1, new_oper, coeff2 = applynoncummuting(gate, operator, theta, coefficient; kwargs...)

    operator_dict[operator] = coeff1
    new_operator_dict[new_oper] = coeff2

    return
end

### Clifford gates

"""
    applygatetoone!(gate::CliffordGate, operator, coefficient, theta, operator_dict, new_operator_dict, args...; kwargs...)

Overload of `applygatetoone!` for `CliffordGate` gates.
Simplified logic for readability.
"""
@inline function applygatetoone!(gate::CliffordGate, operator, coefficient, theta, operator_dict, new_operator_dict, args...; kwargs...)

    op, coeff = apply(gate, operator, theta, coefficient; kwargs...)
    new_operator_dict[op] = coeff

    return
end

### Amplitude Damping Noise
"""
    applygatetoall!(gate::PauliGateUnion, thetas, operator_dict, new_operator_dict, args...; kwargs...)

Overload of `applygatetoall!` for `AmplitudeDampingNoise` gates.
Both `operator_dict` and `new_operator_dict` contain operators which will be merged later.
"""
function applygatetoall!(gate::AmplitudeDampingNoise, theta, operator_dict, new_operator_dict, args...; kwargs...)
    # TODO: there is a lot of code duplication. Can we write a more general function? 

    for (operator, coeff) in operator_dict
        applygatetoone!(gate, operator, coeff, theta, operator_dict, new_operator_dict; kwargs...)
    end

    return operator_dict, new_operator_dict
end

"""
    applygatetoone!(gate::PauliGateUnion, operator, coefficient, theta, operator_dict, new_operator_dict, args...; kwargs...)

Overload of `applygatetoone!` for `AmplitudeDampingNoise` gates.
Checks for whether `gate` will cause splitting and has tailored logic.
"""
@inline function applygatetoone!(gate::AmplitudeDampingNoise, operator, coefficient, theta, operator_dict, new_operator_dict, args...; kwargs...)

    if actsdiagonally(gate, operator)
        operator, coeff = diagonalapply(gate, operator, theta, coefficient; kwargs...)
        operator_dict[operator] = coeff
        return
    end

    operator, coeff1, new_oper, coeff2 = splitapply(gate, operator, theta, coefficient; kwargs...)

    operator_dict[operator] = coeff1
    new_operator_dict[new_oper] = coeff2

    return
end

### Frozen Gates
"""
    applygatetoall!(gate::FrozenGate, thetas, operator_dict, new_operator_dict, args...; kwargs...)

Overload of `applygatetoall!` for `FrozenGate`s. Re-directs to `applygatetoall!` for the wrapped `FrozenGate.gate`.
"""
function applygatetoall!(gate::FrozenGate, theta, operator_dict, new_operator_dict, args...; kwargs...)
    return applygatetoall!(gate.gate, gate.parameter, operator_dict, new_operator_dict, args...; kwargs...)
end