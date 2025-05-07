### paulifreqtracker.jl
##
# PauliFreqTracker type and methods.
# It records the behavior at PauliRotation gates, i.e., the number of times it received a sin or cos factor, and the total number of branchings/splits.
# These path properties can be used for truncations.
# By default, we support `max_freq` and `max_nsins` truncations if the coefficients are of type `PauliFreqTracker`.
##
###

"""
    PauliFreqTracker(coeff::Number, nsins::Int, ncos::Int, freq::Int)

Wrapper type for numerical coefficients in Pauli propagation that records 
the number of sin and cos factors applied via a `PauliRotation` gate, and the so-called frequency, which is their sum.
It appears redundant but these three properties need to be tracked separately because of how merging affects them.
"""
struct PauliFreqTracker{T<:Number} <: PathProperties
    coeff::T
    nsins::Int
    ncos::Int
    freq::Int
end

"""
    PauliFreqTracker(coeff::Number)

Constructor for `PauliFreqTracker` from only a coefficient.
Initializes `nsins`, `ncos`, and `freq` to zero.
"""
PauliFreqTracker(coeff::Number) = PauliFreqTracker(float(coeff), 0, 0, 0)

### Specializations for PauliRotations that incremet the nsins, ncos, and freq

# Overload of `applytoall!` for `PauliRotation` gates acting onto Pauli sums with `PathProperties` coefficients. 
function applytoall!(gate::PauliRotation, theta, psum::PauliSum{TT,PProp}, aux_psum; kwargs...) where {TT<:PauliStringType,PProp<:PathProperties}
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

## Specializations for PauliRotations that increment the nsins, ncos, and freq
# can be used by all PathProperties types that have the necessary fields `ncos`, `nsins`, and `freq`
function splitapply(gate::MaskedPauliRotation, pstr::PauliStringType, coeff::PProp, theta; kwargs...) where {PProp<:PathProperties}
    # increments ncos and freq field if applicable
    coeff1 = _applycos(coeff, theta; kwargs...)
    new_pstr, sign = getnewpaulistring(gate, pstr)
    # increments nsins and freq field if applicable
    coeff2 = _applysin(coeff, theta, sign; kwargs...)

    return pstr, coeff1, new_pstr, coeff2
end

# These also work for other PathProperties types that have a `coeff` field defined
# Multiply sin(theta) * sign to the `coeff` field of a `PathProperties` object.
# Increments the `nsins` and `freq` fields by 1 if applicable.
function _applysin(pth::PProp, theta, sign=1; kwargs...) where {PProp<:PathProperties}
    fields = fieldnames(PProp)

    if :coeff ∉ fields
        throw(
            "The $(PProp) object does not have a field `coeff` to use the `_applysin` operation. " *
            "Consider defining _applysin(pth::$(PProp), theta, sign; kwargs...)"
        )
    end

    function updateval(val, field)
        if field == :coeff
            # apply sin to the `coeff` field
            return val * sin(theta) * sign
        elseif field == :nsins
            # increment the `nsins` field
            return val + 1
        elseif field == :freq
            # increment the `freq` field
            return val + 1
        else
            return val
        end
    end

    return PProp((updateval(getfield(pth, field), field) for field in fields)...)
end


# Multiply cos(theta) * sign to the `coeff` field of a `PathProperties` object.
# Increments the `ncos` and `freq` fields by 1 if applicable.
function _applycos(pth::PProp, theta, sign=1; kwargs...) where {PProp<:PathProperties}
    fields = fieldnames(PProp)

    if :coeff ∉ fields
        throw(
            "The $(PProp) object does not have a field `coeff` to use the `_applysin` operation. " *
            "Consider defining _applycos(pth::$(PProp), theta, sign; kwargs...)"
        )
    end

    function updateval(val, field)
        if field == :coeff
            # apply cos to the `coeff` field
            return val * cos(theta) * sign
        elseif field == :ncos
            # increment the `ncos` field
            return val + 1
        elseif field == :freq
            # increment the `freq` field
            return val + 1
        else
            return val
        end
    end

    return PProp((updateval(getfield(pth, field), field) for field in fields)...)
end