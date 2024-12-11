## This file contains specialized functions for some of our gates. 

### PAULI GATES
"""
    applygatetoone!(gate::PauliRotationUnion, pstr, coefficient, theta, psum, aux_psum, args...; kwargs...)

Overload of `applygatetoone!` for `PauliRotation` and `FastPauliRotation` gates. 
Checks for commutation of `gate` and `pstr`, and applies the gate to the Pauli string if they don't.
"""
@inline function applygatetoone!(gate::PauliRotationUnion, pstr, coefficient, theta, psum, aux_psum, args...; kwargs...)

    if commutes(gate, pstr)
        return
    end

    pstr, coeff1, new_pstr, coeff2 = applynoncummuting(gate, pstr, theta, coefficient; kwargs...)

    psum[pstr] = coeff1
    aux_psum[new_pstr] = coeff2

    return
end

### Clifford gates

"""
    applygatetoall!(gate::CliffordGate, theta, psum, aux_psum, args...; kwargs...)

Overload of `applygatetoall!` for `CliffordGate` gates.


"""
function applygatetoall!(gate::CliffordGate, theta, psum, aux_psum, args...; kwargs...)

    for (pstr, coeff) in psum
        applygatetoone!(gate, pstr, coeff, theta, psum, aux_psum; kwargs...)
    end

    # Empty psum because everything was moved into aux_psum. They will later be merged.
    empty!(psum)

    return psum, aux_psum
end

"""
    applygatetoone!(gate::CliffordGate, pstr, coefficient, theta, psum, aux_psum, args...; kwargs...)

Overload of `applygatetoone!` for `CliffordGate` gates.
Does not delete from `psum`, here and instead empties it in  `applygatetoall!`
"""
@inline function applygatetoone!(gate::CliffordGate, pstr, coefficient, theta, psum, aux_psum, args...; kwargs...)

    new_pstr, coeff = apply(gate, pstr, theta, coefficient; kwargs...)
    aux_psum[new_pstr] = coeff

    return
end

### Amplitude Damping Noise
"""
    applygatetoone!(gate::AmplitudeDampingNoise, pstr, coefficient, theta, psum, aux_psum, args...; kwargs...)

Overload of `applygatetoone!` for `AmplitudeDampingNoise` gates.
Checks for whether `gate` will cause splitting and has tailored logic.
"""
@inline function applygatetoone!(gate::AmplitudeDampingNoise, pstr, coefficient, theta, psum, aux_psum, args...; kwargs...)

    if actsdiagonally(gate, pstr)
        pstr, coeff = diagonalapply(gate, pstr, theta, coefficient; kwargs...)
        psum[pstr] = coeff
        return
    end

    pstr, coeff1, new_pstr, coeff2 = splitapply(gate, pstr, theta, coefficient; kwargs...)

    psum[pstr] = coeff1
    aux_psum[new_pstr] = coeff2

    return
end

### Frozen Gates
"""
    applygatetoall!(gate::FrozenGate, thetas, psum, aux_psum, args...; kwargs...)

Overload of `applygatetoall!` for `FrozenGate`s. Re-directs to `applygatetoall!` for the wrapped `FrozenGate.gate`.
"""
function applygatetoall!(gate::FrozenGate, theta, psum, aux_psum, args...; kwargs...)
    return applygatetoall!(gate.gate, gate.parameter, psum, aux_psum, args...; kwargs...)
end