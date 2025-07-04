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
    applytoall!(gate::PauliRotation, theta, psum, aux_psum; kwargs...)

Overload of `applytoall!` for `PauliRotation` gates. 
It fixes the type-instability of the `apply()` function and reduces moving Pauli strings between `psum` and `aux_psum`.
`psum` and `aux_psum` are merged later.
"""
function applytoall!(gate::PauliRotation, theta, psum, aux_psum; kwargs...)
    # turn the (potentially) PauliRotation gate into a MaskedPauliRotation gate
    # this allows for faster operations
    gate = _tomaskedpaulirotation(gate, paulitype(psum))
    applytoall!(gate, theta, psum, aux_psum; kwargs...)
end

function applytoall!(gate::MaskedPauliRotation, theta, psum, aux_psum; kwargs...)
    # pre-compute the sine and cosine values because the are used for every Pauli string that does not commute with the gate
    cos_val = cos(theta)
    sin_val = sin(theta)
    # loop over all Pauli strings and their coefficients in the Pauli sum
    for (pstr, coeff) in psum

        if commutes(gate, pstr)
            # if the gate commutes with the pauli string, do nothing
            continue
        end

        # else we know the gate will split th Pauli string into two
        coeff1 = coeff * cos_val
        new_pstr, sign = getnewpaulistring(gate, pstr)
        coeff2 = coeff * sin_val * sign

        # set the coefficient of the original Pauli string
        set!(psum, pstr, coeff1)

        # set the coefficient of the new Pauli string in the aux_psum
        # we can set the coefficient because PauliRotations create non-overlapping new Pauli strings
        set!(aux_psum, new_pstr, coeff2)
    end

    return
end


function getnewpaulistring(gate::MaskedPauliRotation, pstr::PauliStringType)
    new_pstr = _bitpaulimultiply(gate.generator_mask, pstr)

    # this counts the exponent of the imaginary unit in the new Pauli string
    im_count = _calculatesignexponent(gate.generator_mask, pstr)

    # now, instead of computing im^im_count followed by another im factor from the gate rules,
    # we do this in one step via a cheeky trick:
    sign = (im_count & 2) - 1
    # this is equivalent to sign = real( im * im^im_count)

    return new_pstr, sign
end


### Clifford gates
"""
    applyandadd!(gate::CliffordGate, pstr, coeff, theta, output_psum; kwargs...)

Overload of `applyandadd!` for `CliffordGate` gates.
Use `set!()` instead of `add!()` because Clifford gates create non-overlapping Pauli strings.
`applytoall!` does not need to be adapted.
"""
@inline function applyandadd!(gate::CliffordGate, pstr, coeff, theta, output_psum; kwargs...)

    # TODO: test whether it is significantly faster to get the map_array in applytoall! and pass it here
    new_pstr, new_coeff = apply(gate, pstr, coeff; kwargs...)
    # we can set the coefficient because Cliffords create non-overlapping Pauli strings
    set!(output_psum, new_pstr, new_coeff)

    return
end

"""
    apply(gate::CliffordGate, pstr::PauliStringType, coeff)

Apply a `CliffordGate` to an integer Pauli string and its coefficient. 
"""
function apply(gate::CliffordGate, pstr::PauliStringType, coeff; kwargs...)
    # this array carries the new Paulis + sign for every occuring old Pauli combination
    map_array = clifford_map[gate.symbol]

    qinds = gate.qinds

    # this integer carries the active Paulis on its bits
    lookup_int = getpauli(pstr, qinds)

    # this integer can be used to index into the array returning the new Paulis
    # +1 because Julia is 1-indexed and lookup_int is 0-indexed
    partial_pstr, sign = map_array[lookup_int+1]

    # insert the bits of the new Pauli into the old Pauli
    pstr = setpauli(pstr, partial_pstr, qinds)

    coeff *= sign

    # usually return format would be ((pstr, coeff), ), but we also overload applyandadd!() to accept this
    return pstr, coeff
end

### Pauli Noise
"""
    applytoall!(gate::PauliNoise, p, psum, aux_psum; kwargs...)

Overload of `applytoall!` for `PauliNoise` gates with noise strength `p`. 
It changes the coefficients in-place and does not require the `aux_psum`, which stays empty.
"""
function applytoall!(gate::PauliNoise, p, psum, aux_psum; kwargs...)
    # check that the noise strength is in the correct range
    _check_noise_strength(PauliNoise, p)

    # loop over all Pauli strings and their coefficients in the Pauli sum
    for (pstr, coeff) in psum

        # the Pauli on the site that the noise acts on
        pauli = getpauli(pstr, gate.qind)

        # `isdamped` is defined in noisechannels.jl for each Pauli noise channel
        # I Paulis are never damped, but the others vary
        if !isdamped(gate, pauli)
            continue
        end

        new_coeff = coeff * (1 - p)
        # change the coefficient in psum, don't move anything to aux_psum
        set!(psum, pstr, new_coeff)
    end

    return
end

### Amplitude Damping Noise
"""
    applytoall!(gate::AmplitudeDampingNoise, theta, psum, aux_psum; kwargs...)

Overload of `applytoall!` for `AmplitudeDampingNoise` gates. 
It fixes the type-instability of the apply() function and reduces moving Pauli strings between psum and aux_psum.
`psum` and `aux_psum` are merged later.
"""
function applytoall!(gate::AmplitudeDampingNoise, gamma, psum, aux_psum; kwargs...)
    # check that the noise strength is in the correct range
    _check_noise_strength(AmplitudeDampingNoise, gamma)

    # loop over all Pauli strings and their coefficients in the Pauli sum
    for (pstr, coeff) in psum
        pauli = getpauli(pstr, gate.qind)
        if pauli == 0
            # Pauli is I, so the gate does not do anything
            continue

        elseif pauli == 1 || pauli == 2
            # Pauli is X or Y, so the gate will give a sqrt(1-gamma) prefactor
            new_coeff = sqrt(1 - gamma)
            # set the coefficient of the Pauli string in the psum to the new coefficient
            set!(psum, pstr, new_coeff)

        else
            # Pauli is Z, so the gate will split the Pauli string 

            # else we know the gate will split th Pauli string into two
            new_pstr = setpauli(pstr, 0, gate.qind)
            coeff1 = (1 - gamma) * coeff
            coeff2 = gamma * coeff

            # set the coefficient of the original Pauli string
            set!(psum, pstr, coeff1)

            # add the coefficient of the new Pauli string in the aux_psum
            add!(aux_psum, new_pstr, coeff2)

        end
    end

    return
end

## T Gate
"""
    applytoall!(gate::TGate, thetas, psum, aux_psum; kwargs...)

Overload of `applytoall!()` for `TGate(qind)`.
Redirects to a `PauliRotation(:Z, qind)` with angle π/4.
"""
function applytoall!(gate::TGate, thetas, psum, aux_psum; kwargs...)
    return applytoall!(PauliRotation(:Z, gate.qind), π / 4, psum, aux_psum; kwargs...)
end

## TransferMapGate
"""
    apply(gate::TransferMapGate, pstr, coeff)

Apply a `TransferMapGate` to an integer Pauli string and its coefficient.
The outcomes are determined by the `transfer_map` of the gate.
"""
function PauliPropagation.apply(gate::TransferMapGate, pstr, coeff)
    # the Paulis packed into the integer are used to index into the transfer map
    pauli_int = getpauli(pstr, gate.qinds)
    pstrs_and_factors = gate.transfer_map[pauli_int+1]
    # the new pstrs are the new Paulis that need to be set and the coefficients need to be multiplied with the factors
    return Tuple((setpauli(pstr, new_pstr, gate.qinds), coeff * factor) for (new_pstr, factor) in pstrs_and_factors)
end

### Frozen Gates
"""
    applytoall!(gate::FrozenGate, thetas, psum, aux_psum; kwargs...)

Overload of `applytoall!` for `FrozenGate`s. Re-directs to `applytoall!` for the wrapped `FrozenGate.gate` with the frozen parameter.
"""
function applytoall!(gate::FrozenGate, theta, psum, aux_psum; kwargs...)
    return applytoall!(gate.gate, gate.parameter, psum, aux_psum; kwargs...)
end
