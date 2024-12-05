## This file contains specialized functions for some of our gates. 

### PAULI GATES
"""
    applygatetoall!(gate::PauliGateUnion, thetas, psum, second_psum, args...; kwargs...)

Overload of `applygatetoall!` for `PauliGate` and `FastPauliGate` gates.
Both `psum` and `second_psum` contain Pauli strings which will be merged later.
"""
function applygatetoall!(gate::PauliGateUnion, theta, psum, second_psum, args...; kwargs...)
    # TODO: there is a lot of code duplication. Can we write a more general function?

    for (pstr, coeff) in psum
        applygatetoone!(gate, pstr, coeff, theta, psum, second_psum; kwargs...)
    end

    return psum, second_psum  # don't swap psums around
end

"""
    applygatetoone!(gate::PauliGateUnion, pstr, coefficient, theta, psum, second_psum, args...; kwargs...)

Overload of `applygatetoone!` for `PauliGate` and `FastPauliGate` gates. 
Checks for commutation of `gate` and `pstr`, and applies the gate to the Pauli string if they don't.
"""
@inline function applygatetoone!(gate::PauliGateUnion, pstr, coefficient, theta, psum, second_psum, args...; kwargs...)

    if commutes(gate, pstr)
        return
    end

    pstr, coeff1, new_pstr, coeff2 = applynoncummuting(gate, pstr, theta, coefficient; kwargs...)

    psum[pstr] = coeff1
    second_psum[new_pstr] = coeff2

    return
end

### Clifford gates

"""
    applygatetoone!(gate::CliffordGate, pstr, coefficient, theta, psum, second_psum, args...; kwargs...)

Overload of `applygatetoone!` for `CliffordGate` gates.
Simplified logic for readability.
"""
@inline function applygatetoone!(gate::CliffordGate, pstr, coefficient, theta, psum, second_psum, args...; kwargs...)

    new_pstr, coeff = apply(gate, pstr, theta, coefficient; kwargs...)
    second_psum[new_pstr] = coeff

    return
end

### Amplitude Damping Noise
"""
    applygatetoall!(gate::AmplitudeDampingNoise, thetas, psum, second_psum, args...; kwargs...)

Overload of `applygatetoall!` for `AmplitudeDampingNoise` gates.
Both `psum` and `second_psum` contain Pauli strings which will be merged later.
"""
function applygatetoall!(gate::AmplitudeDampingNoise, theta, psum, second_psum, args...; kwargs...)
    # TODO: there is a lot of code duplication. Can we write a more general function? 

    for (pstr, coeff) in psum
        applygatetoone!(gate, pstr, coeff, theta, psum, second_psum; kwargs...)
    end

    return psum, second_psum
end

"""
    applygatetoone!(gate::AmplitudeDampingNoise, pstr, coefficient, theta, psum, second_psum, args...; kwargs...)

Overload of `applygatetoone!` for `AmplitudeDampingNoise` gates.
Checks for whether `gate` will cause splitting and has tailored logic.
"""
@inline function applygatetoone!(gate::AmplitudeDampingNoise, pstr, coefficient, theta, psum, second_psum, args...; kwargs...)

    if actsdiagonally(gate, pstr)
        pstr, coeff = diagonalapply(gate, pstr, theta, coefficient; kwargs...)
        psum[pstr] = coeff
        return
    end

    pstr, coeff1, new_pstr, coeff2 = splitapply(gate, pstr, theta, coefficient; kwargs...)

    psum[pstr] = coeff1
    second_psum[new_pstr] = coeff2

    return
end

### Frozen Gates
"""
    applygatetoall!(gate::FrozenGate, thetas, psum, second_psum, args...; kwargs...)

Overload of `applygatetoall!` for `FrozenGate`s. Re-directs to `applygatetoall!` for the wrapped `FrozenGate.gate`.
"""
function applygatetoall!(gate::FrozenGate, theta, psum, second_psum, args...; kwargs...)
    return applygatetoall!(gate.gate, gate.parameter, psum, second_psum, args...; kwargs...)
end