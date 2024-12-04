"""
Abstract type for wrapping coefficients and record custom path properties
"""
abstract type PathProperties end

"""
Pretty print for PathProperties
"""
Base.show(io::IO, pth::PathProperties) = print(io, "PathProperties($(typeof(pth.coeff)), nsins=$(pth.nsins), ncos=$(pth.ncos), freq=$(pth.freq))")

import Base: copy
"""
Copy a `PathProperties` object. Does not deepcopy, i.e., the fields are not copied.
"""
function copy(pth::T) where {T<:PathProperties}
    fields = fieldnames(T)
    return T((getfield(pth, field) for field in fields)...)
end


import Base: *
"""
Multiplication of the `coeff` field in a `PathProperties` object with a number.
Requires that the `PathProperties` object has a `coeff` field as the first field.
"""
function *(pth::T, val::Number) where {T<:PathProperties}
    fields = fieldnames(T)
    return T(getfield(pth, :coeff) * val, (getfield(pth, fields[ii]) for ii in 2:length(fields))...)
end

"""
Multiplication of a `PathProperties` object with a number.
Requires that the `PathProperties` object has a `coeff` field as the first field which will be multiplied.
"""
function *(val::Number, pth::T) where {T<:PathProperties}
    return pth * val
end

import Base: +
"""
Addition of two `PathProperties` objects of equal concrete type.
Adds the `coeff` fields and takes the minimum of the other fields.
Requires that the `PathProperties` object has a `coeff` field as the first field.
"""
function +(pth1::T, pth2::T) where {T<:PathProperties}
    fields = fieldnames(T)
    return T(getfield(pth1, :coeff) + getfield(pth2, :coeff), (min(getfield(pth1, fields[ii]), getfield(pth2, fields[ii])) for ii in 2:length(fields))...)
end

"""
    NumericPathProperties(coeff::Number, nsins::Int, ncos::Int, freq::Int)

Wrapper type for numerical coefficients in Pauli propagation that records 
the number of sin and cos terms applied, and the so-called frequency, which is their sum.
It appears redundant but these three properties need to be tracked separately because of how merging affects them.
"""
struct NumericPathProperties{T<:Number} <: PathProperties
    coeff::T
    nsins::Int
    ncos::Int
    freq::Int
end

"""
    NumericPathProperties(coeff::Number)

Constructor for `NumericPathProperties` from only a coefficient.
Initializes `nsins`, `ncos`, and `freq` to zero.
"""
NumericPathProperties(coeff::Number) = NumericPathProperties(float(coeff), 0, 0, 0)

# TODO: More general show method for general fields.
"""
Pretty print for NumericPathProperties
"""
Base.show(io::IO, pth::NumericPathProperties) = print(io, "NumericPathProperties($(pth.coeff), nsins=$(pth.nsins), ncos=$(pth.ncos), freq=$(pth.freq))")


"""
    numcoefftype(pth::PathProperties)

Return the type of the coefficient in a `PathProperties` object.
"""
function numcoefftype(pth::PathProperties)
    return typeof(pth.coeff)
end

"""
    getnumcoeff(val::PathProperties)

Get the numerical coefficient of a `PathProperties` wrapper.
"""
function getnumcoeff(val::PathProperties)
    return val.coeff
end

"""
    wrapcoefficients(pstr::PauliString)

Wrap the coefficient of a `PauliString` into `NumericPathProperties`.
"""
function wrapcoefficients(pstr::PauliString)
    # by default wrap into `NumericPathProperties`
    return wrapcoefficients(pstr, NumericPathProperties)
end

"""
    wrapcoefficients(pstr::PauliString, PathPropertiesType::Type{PP}) where {PP<:PathProperties}

Wrap the coefficient of a `PauliString` into a custom `PathProperties` type.
For anything that is not natively supported by the library, you can subtype `PathProperties`.
A one-argument constructor of the custom `PathProperties` type from a coefficient must be defined.
"""
function wrapcoefficients(pstr::PauliString, PathPropertiesType::Type{PP}) where {PP<:PathProperties}
    # the one-argument constructor of your PathProperties type must be defined
    return PauliString(pstr.nqubits, term(pstr), PathPropertiesType(pstr.coeff))
end

"""
    wrapcoefficients(psum::PauliSum)

Wrap the coefficients of a `PauliSum` into `NumericPathProperties`.
"""
function wrapcoefficients(psum::PauliSum)
    # by default wrap into `NumericPathProperties`
    return wrapcoefficients(psum, NumericPathProperties)
end

# TODO: This is not type stable
"""
    wrapcoefficients(psum::PauliSum, PathPropertiesType::Type{PP}) where {PP<:PathProperties}

Wrap the coefficients of a `PauliSum` into a custom `PathProperties` type.
For anything that is not natively supported by the library, you can subtype `PathProperties`.
A one-argument constructor of the custom `PathProperties` type from a coefficient must be defined.
"""
function wrapcoefficients(psum::PauliSum, PathPropertiesType::Type{PP}) where {PP<:PathProperties}
    return PauliSum(psum.nqubits, Dict(op => PathPropertiesType(coeff) for (op, coeff) in psum.op_dict))
end