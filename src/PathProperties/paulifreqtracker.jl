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
"""
    splitapply(gate::MaskedPauliRotation, pstr::PauliStringType, coeff, theta; kwargs...)

Apply a `MaskedPauliRotation` with an angle `theta` and a coefficient `coeff` to an integer Pauli string,
assuming that the gate does not commute with the Pauli string.
Returns two pairs of (pstr, coeff) as one tuple.
Currently `kwargs` are passed to `applycos` and `applysin` for the Surrogate.
"""
function splitapply(gate::MaskedPauliRotation, pstr::PauliStringType, coeff::PauliFreqTracker, theta; kwargs...)
    coeff1 = _applycos(coeff, theta; kwargs...)
    new_pstr, sign = getnewpaulistring(gate, pstr)
    coeff2 = _applysin(coeff, theta, sign; kwargs...)

    return pstr, coeff1, new_pstr, coeff2
end

# These also work for other PathProperties types that use these function
"""
    _applysin(pth::PathProperties, theta; sign=1, kwargs...)

Multiply sin(theta) * sign to the `coeff` field of a `PathProperties` object.
Increments the `nsins` and `freq` fields by 1 if applicable.
"""
function _applysin(pth::PProp, theta, sign=1; kwargs...) where {PProp<:PathProperties}
    fields = fieldnames(PProp)

    function updateval(val, field)
        if field == :coeff
            # apply sin to the `coeff` field
            # TODO: Should we even recursively call _applysin and _applycos?
            return _applysin(val, theta, sign; kwargs...)
        elseif field == :nsins
            return val + 1
        elseif field == :freq
            return val + 1
        else
            return val
        end
    end
    return PProp((updateval(getfield(pth, field), field) for field in fields)...)
end

"""
    _applycos(pth::PathProperties, theta; sign=1, kwargs...)

Multiply cos(theta) * sign to the `coeff` field of a `PathProperties` object.
Increments the `ncos` and `freq` fields by 1 if applicable.
"""
function _applycos(pth::PProp, theta, sign=1; kwargs...) where {PProp<:PathProperties}
    fields = fieldnames(PProp)

    function updateval(val, field)
        if field == :coeff
            # apply cos to the `coeff` field
            return _applycos(val, theta, sign; kwargs...)
        elseif field == :ncos
            return val + 1
        elseif field == :freq
            return val + 1
        else
            return val
        end
    end
    return PProp((updateval(getfield(pth, field), field) for field in fields)...)
end

_applycos(val::Number, theta, sign=1; kwargs...) = val * cos(theta) * sign
_applysin(val::Number, theta, sign=1; kwargs...) = val * sin(theta) * sign
