###
##
# This file contains specialized functions for some of our gates.
# We overload `applyandadd!()` to fix potential type-instabilities in `apply()` if the number of returned Pauli strings is not fixed.
# We overload `applytoall!()` to reduce unnecessarily moving Pauli strings between `psum` and `aux_psum`.
# This usually also fixes potential type-instabilities in `apply()`.
# Both functions can be overloaded if needed.
##
###

### PAULI GATES
"""
    applytoall!(gate::PauliRotationUnion, theta, psum, aux_psum, args...; kwargs...)

Overload of `applytoall!` for `PauliRotation` and `MaskedPauliRotation` gates. 
It fixes the type-instability of the `apply()` function and reduces moving Pauli strings between `psum` and `aux_psum`.
`psum` and `aux_psum` are merged later.
"""
function applytoall!(gate::PauliRotationUnion, theta, psum, aux_psum, args...; kwargs...)
    # turn the (potentially) PauliRotation gate into a MaskedPauliRotation gate
    # this allows for faster operations
    gate = _tomaskedpaulirotation(gate, paulitype(psum))

    # loop over all Pauli strings and their coefficients in the Pauli sum
    for (pstr, coeff) in psum

        if commutes(gate, pstr)
            # if the gate commutes with the pauli string, do nothing
            continue
        end

        # else we know the gate will split th Pauli string into two
        pstr, coeff1, new_pstr, coeff2 = splitapply(gate, pstr, theta, coeff; kwargs...)

        # set the coefficient of the original Pauli string
        set!(psum, pstr, coeff1)

        # set the coefficient of the new Pauli string in the aux_psum
        # we can set the coefficient because PauliRotations create non-overlapping new Pauli strings
        set!(aux_psum, new_pstr, coeff2)
    end

    return
end


### Clifford gates
"""
    applyandadd!(gate::CliffordGate, pstr, coeff, theta, output_psum, args...; kwargs...)

Overload of `applyandadd!` for `CliffordGate` gates.
Use `set!()` instead of `add!()` because Clifford gates create non-overlapping Pauli strings.
`applytoall!` does not need to be adapted.
"""
@inline function applyandadd!(gate::CliffordGate, pstr, coeff, theta, output_psum, args...; kwargs...)

    # TODO: test whether it is significantly faster to get the map_array in applytoall! and pass it here
    new_pstr, new_coeff = apply(gate, pstr, theta, coeff; kwargs...)
    # we can set the coefficient because Cliffords create non-overlapping Pauli strings
    set!(output_psum, new_pstr, new_coeff)

    return
end

"""
    apply(gate::CliffordGate, pstr::PauliStringType, coefficient=1.0)

Apply a `CliffordGate` to an integer Pauli string and an optional coefficient. 
"""
function apply(gate::CliffordGate, pstr::PauliStringType, coefficient=1.0; kwargs...)
    map_array = clifford_map[gate.symbol]
    return applywithmap(gate, pstr, coefficient, map_array)
end

### Pauli Noise
"""
    applytoall!(gate::PauliNoise, theta, psum, aux_psum, args...; kwargs...)

Overload of `applytoall!` for `PauliNoise` gates. 
It changes the coefficients in-place and does not require the `aux_psum`, which stays empty.
"""
function applytoall!(gate::PauliNoise, theta, psum, aux_psum, args...; kwargs...)

    # loop over all Pauli strings and their coefficients in the Pauli sum
    for (pstr, coeff) in psum
        if getpauli(pstr, gate.qind) == 0
            # Pauli is I, so the gate does not do anything
            continue
        end

        # apply the Pauli noise, which will reduce the coefficient
        pstr, new_coeff = apply(gate, pstr, theta, coeff; kwargs...)

        # set the coefficient of the Pauli string in the psum to the new coefficient
        set!(psum, pstr, new_coeff)
    end

    return
end


"""
    apply(gate::DepolarizingNoise, pstr::PauliStringType, p, coefficient=1.0)

Apply a depolarizing noise channel to an integer Pauli string `pstr` with noise strength `p`.
Physically `p` is restricted to the range `[0, 1]`.
A coefficient of the Pauli string can optionally be passed as `coefficient`.
"""
function apply(gate::PauliNoise, pstr::PauliStringType, p, coefficient; kwargs...)
    pauli = getpauli(pstr, gate.qind)

    if isdamped(gate, pauli) # this function is defined in the noisechannels.jl file for each Pauli noise channel
        coefficient *= (1 - p)
    end

    return pstr, coefficient
end


### Amplitude Damping Noise
"""
    applytoall!(gate::AmplitudeDampingNoise, theta, psum, aux_psum, args...; kwargs...)

Overload of `applytoall!` for `AmplitudeDampingNoise` gates. 
It fixes the type-instability of the apply() function and reduces moving Pauli strings between psum and aux_psum.
`psum` and `aux_psum` are merged later.
"""
function applytoall!(gate::AmplitudeDampingNoise, theta, psum, aux_psum, args...; kwargs...)

    # loop over all Pauli strings and their coefficients in the Pauli sum
    for (pstr, coeff) in psum
        pauli = getpauli(pstr, gate.qind)
        if pauli == 0
            # Pauli is I, so the gate does not do anything
            continue
        elseif pauli == 1 || pauli == 2
            # Pauli is X or Y, so the gate will give a sqrt(1-gamma) prefactor
            pstr, new_coeff = diagonalapply(gate, pstr, theta, coeff; kwargs...)
            # set the coefficient of the Pauli string in the psum to the new coefficient
            set!(psum, pstr, new_coeff)
        else
            # Pauli is Z, so the gate will split the Pauli string 

            # else we know the gate will split th Pauli string into two
            pstr, coeff1, new_pstr, coeff2 = splitapply(gate, pstr, theta, coeff; kwargs...)

            # set the coefficient of the original Pauli string
            set!(psum, pstr, coeff1)

            # add the coefficient of the new Pauli string in the aux_psum
            add!(aux_psum, new_pstr, coeff2)

        end
    end

    return
end


### Frozen Gates
"""
    applytoall!(gate::FrozenGate, thetas, psum, aux_psum, args...; kwargs...)

Overload of `applytoall!` for `FrozenGate`s. Re-directs to `applytoall!` for the wrapped `FrozenGate.gate` with the frozen parameter.
"""
function applytoall!(gate::FrozenGate, theta, psum, aux_psum, args...; kwargs...)
    return applytoall!(gate.gate, gate.parameter, psum, aux_psum, args...; kwargs...)
end

# TODO: to we need the apply() function at all?
# """
#     apply(frozen_gate::FrozenGate, pstr, theta, args...)

# Apply a `FrozenGate` to a Pauli string `pstr` with the parameter `FrozenGate.parameter`.
# The passed `theta` will be ignored.
# """
# function apply(frozen_gate::FrozenGate, pstr, theta, coefficient=1.0; kwargs...)
#     return apply(frozen_gate.gate, pstr, frozen_gate.parameter, coefficient; kwargs...)
# end