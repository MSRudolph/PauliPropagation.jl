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
        pstr, coeff1, new_pstr, coeff2 = splitapply(gate, pstr, coeff, theta; kwargs...)

        # set the coefficient of the original Pauli string
        set!(psum, pstr, coeff1)

        # set the coefficient of the new Pauli string in the aux_psum
        # we can set the coefficient because PauliRotations create non-overlapping new Pauli strings
        set!(aux_psum, new_pstr, coeff2)
    end

    return
end

"""
    splitapply(gate::PauliRotation, pstr::PauliStringType, coeff, theta; kwargs...)

Apply a `PauliRotation` with an angle `theta` and a coefficient `coeff` to an integer Pauli string,
assuming that the gate does not commute with the Pauli string.
Returns two pairs of (pstr, coeff) as one tuple.
Currently `kwargs` are passed to `applycos` and `applysin` for the Surrogate.
"""
function splitapply(gate::PauliRotationUnion, pstr::PauliStringType, coeff, theta; kwargs...)
    coeff1 = _applycos(coeff, theta; kwargs...)
    new_pstr, sign = getnewpaulistring(gate, pstr)
    coeff2 = _applysin(coeff, theta; sign=sign, kwargs...)

    return pstr, coeff1, new_pstr, coeff2
end

"""
    getnewpaulistring(gate::PauliRotation, pstr::PauliStringType)

Get the new Pauli string after applying a `PauliRotation` to an integer Pauli string,
as well as the corresponding ±1 coefficient.
"""
function getnewpaulistring(gate::PauliRotation, pstr::PauliStringType)
    new_pstr = copy(pstr)

    total_sign = 1  # this coefficient will be imaginary
    for (qind, gate_sym) in zip(gate.qinds, gate.symbols)
        sign, new_partial_str = pauliprod(gate_sym, getpauli(new_pstr, qind))
        total_sign *= sign
        new_pstr = setpauli(new_pstr, new_partial_str, qind)
    end
    return new_pstr, real(1im * total_sign)
end

"""
    getnewpaulistring(gate::MaskedPauliRotation, pstr::PauliStringType)

Get the new Pauli string after applying a `MaskedPauliRotation` to an integer Pauli string,
as well as the corresponding ±1 coefficient.
"""
function getnewpaulistring(gate::MaskedPauliRotation, pstr::PauliStringType)
    sign, new_pstr = pauliprod(gate.generator_mask, pstr, gate.qinds)
    return new_pstr, real(1im * sign)
end

# NOTE: The reason we have the following functions defined is so PathProperties and Surrogate can overload them
"""
    _applysin(old_coeff::Number, theta; sign=1, kwargs...)

Multiply a numerical coefficient with sin(theta) * sign.
"""
function _applysin(old_coeff::Number, theta; sign=1, kwargs...)
    return old_coeff * sin(theta) * sign
end

"""
    _applycos(old_coeff::Number, theta; sign=1, kwargs...)

Multiply a numerical coefficient with cos(theta) * sign.
"""
function _applycos(old_coeff::Number, theta; sign=1, kwargs...)
    return old_coeff * cos(theta) * sign
end

function _incrementcosandfreq(coeff)
    return coeff
end

function _incrementsinandfreq(coeff)
    return coeff
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
    new_pstr, new_coeff = apply(gate, pstr, coeff; kwargs...)
    # we can set the coefficient because Cliffords create non-overlapping Pauli strings
    set!(output_psum, new_pstr, new_coeff)

    return
end

"""
    apply(gate::CliffordGate, pstr::PauliStringType, coeff)

Apply a `CliffordGate` to an integer Pauli string and an optional coefficient. 
"""
function apply(gate::CliffordGate, pstr::PauliStringType, coeff; kwargs...)
    map_array = clifford_map[gate.symbol]
    return applywithmap(gate, pstr, coeff, map_array)
end

"""
    applywithmap(gate::CliffordGate, pstr::PauliStringType, coefficient, map_array)

Apply a `CliffordGate` to an integer Pauli string and a coefficient 
using the a `map_array` corresponding to the `CliffordGate`.
"""
function applywithmap(gate::CliffordGate, pstr::PauliStringType, coeff, map_array; kwargs...)
    qinds = gate.qinds

    lookup_int = _extractlookupop(pstr, qinds)
    sign, partial_pstr = map_array[lookup_int+1]  # +1 because Julia is 1-indexed and lookup_int is 0-indexed
    pstr = _insertnewpaulis!(pstr, partial_pstr, qinds)

    coeff = _multiplysign(coeff, sign)
    return pstr, coeff
end

function _extractlookupop(lookup_int::PauliStringType, qinds)
    partial_pstr = typeof(lookup_int)(0)
    for ii in eachindex(qinds)
        partial_pstr = setpauli(partial_pstr, getpauli(lookup_int, qinds[ii]), ii)
    end
    return partial_pstr
end

function _insertnewpaulis!(pstr::PauliStringType, partial_pstr::PauliStringType, qinds)
    for ii in eachindex(qinds)
        pstr = setpauli(pstr, getpauli(partial_pstr, ii), qinds[ii])
    end
    return pstr
end

function _multiplysign(coefficient, sign)
    return coefficient * sign
end

### Pauli Noise
"""
    applytoall!(gate::PauliNoise, p, psum, aux_psum, args...; kwargs...)

Overload of `applytoall!` for `PauliNoise` gates with noise strength `p`. 
It changes the coefficients in-place and does not require the `aux_psum`, which stays empty.
"""
function applytoall!(gate::PauliNoise, p, psum, aux_psum, args...; kwargs...)

    # loop over all Pauli strings and their coefficients in the Pauli sum
    for (pstr, coeff) in psum
        if getpauli(pstr, gate.qind) == 0
            # Pauli is I, so the gate does not do anything
            continue
        end

        # apply the Pauli noise, which will reduce the coefficient
        pstr, new_coeff = apply(gate, pstr, coeff, p; kwargs...)

        # set the coefficient of the Pauli string in the psum to the new coefficient
        set!(psum, pstr, new_coeff)
    end

    return
end


"""
    apply(gate::PauliNoise, pstr::PauliStringType, coeff, p)

Apply a `PauliNoise` channel to an integer Pauli string `pstr` with noise strength `p`.
Physically `p` is restricted to the range `[0, 1]` and updates the coefficient `coeff` via `coeff*(1-p)` if damped.
"""
function apply(gate::PauliNoise, pstr::PauliStringType, coeff, p; kwargs...)
    pauli = getpauli(pstr, gate.qind)

    if isdamped(gate, pauli) # this function is defined in the noisechannels.jl file for each Pauli noise channel
        coeff *= (1 - p)
    end

    return pstr, coeff
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
            pstr, new_coeff = diagonalapply(gate, pstr, coeff, theta; kwargs...)
            # set the coefficient of the Pauli string in the psum to the new coefficient
            set!(psum, pstr, new_coeff)
        else
            # Pauli is Z, so the gate will split the Pauli string 

            # else we know the gate will split th Pauli string into two
            pstr, coeff1, new_pstr, coeff2 = splitapply(gate, pstr, coeff, theta; kwargs...)

            # set the coefficient of the original Pauli string
            set!(psum, pstr, coeff1)

            # add the coefficient of the new Pauli string in the aux_psum
            add!(aux_psum, new_pstr, coeff2)

        end
    end

    return
end

"""
    actsdiagonally(gate::AmplitudeDampingNoise, pstr::PauliStringType)

Check if the amplitude damping noise channel acts diagonally on the Pauli string `pstr`.
This implies no splitting (acting diagonally) which happens when acting on I, X, and Y.
"""
function actsdiagonally(gate::AmplitudeDampingNoise, pstr::PauliStringType; kwargs...)
    return getpauli(pstr, gate.qind) != 3
end

"""
    diagonalapply(gate::AmplitudeDampingNoise, pstr::PauliStringType, coeff, gamma)

Apply an amplitude damping noise channel to an integer Pauli string `pstr` with noise strength `gamma`.
This is under the assumption that it has been checked that the noise channel acts diagonally on the Pauli string.
Returns a tuple of Pauli string and coefficient.
Physically `gamma` is restricted to the range `[0, 1]`.
A coefficient of the Pauli string can optionally be passed as `coefficient`.
"""
function diagonalapply(gate::AmplitudeDampingNoise, pstr::PauliStringType, coeff, gamma; kwargs...)

    local_pauli = getpauli(pstr, gate.qind)

    if local_pauli != 0  # non-identity Pauli
        coeff *= sqrt(1 - gamma)
    end

    return pstr, coeff
end

"""
    splitapply(gate::AmplitudeDampingNoise, pstr::PauliStringType, coeff, gamma)

Apply an amplitude damping noise channel to an integer Pauli string `pstr` with noise strength `gamma`.
This is under the assumption that it has been checked that the noise channel acts on a Z Pauli and splits.
Returns a tuple of two pairs of Pauli strings and coefficients.
Physically `gamma` is restricted to the range `[0, 1]`.
A coefficient of the Pauli string can optionally be passed as `coefficient`.
"""
function splitapply(gate::AmplitudeDampingNoise, pstr::PauliStringType, coeff, gamma; kwargs...)
    new_pstr = setpauli(pstr, 0, gate.qind)
    return pstr, (1 - gamma) * coeff, new_pstr, gamma * coeff
end

### Frozen Gates
"""
    applytoall!(gate::FrozenGate, thetas, psum, aux_psum, args...; kwargs...)

Overload of `applytoall!` for `FrozenGate`s. Re-directs to `applytoall!` for the wrapped `FrozenGate.gate` with the frozen parameter.
"""
function applytoall!(gate::FrozenGate, theta, psum, aux_psum, args...; kwargs...)
    return applytoall!(gate.gate, gate.parameter, psum, aux_psum, args...; kwargs...)
end