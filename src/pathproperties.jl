"""
Abstract type for wrapping coefficients and record custom path properties
"""
abstract type PathProperties end

"""
Pretty print for PathProperties
"""
Base.show(io::IO, pth::PathProperties) = print(io, "PathProperties($(typeof(pth.coeff)), nsins=$(pth.nsins), ncos=$(pth.ncos), freq=$(pth.freq))")

import Base: *
"""
Multiplication of `PathProperties` with a number modifies coefficient in-place.
"""
function *(pth::PathProperties, val::Number)
    pth.coeff *= val
    return pth
end
import Base: copy
"""
Copy a `PathProperties` object. Does not copy the coefficient `coeff`.
"""
function copy(path_properties::PathProperties)
    return typeof(path_properties)(path_properties.coeff, path_properties.nsins, path_properties.ncos, path_properties.freq)
end

"""
    NumericPathProperties(coeff::Float64, nsins::Int, ncos::Int, freq::Int)

Wrapper type for numerical coefficients in Pauli propagation that records 
the number of sin and cos terms applied, and the so-called frequency, which is their sum.
It appears redundant but these three properties need to be tracked separately because of how merging affects them.
"""
mutable struct NumericPathProperties <: PathProperties
    coeff::Float64
    nsins::Int
    ncos::Int
    freq::Int
end

#TODO: Adapt to non-Float64 coefficients.
"""
    NumericPathProperties(coeff::Number)

Constructor for `NumericPathProperties` from only a coefficient.
Initializes `nsins`, `ncos`, and `freq` to zero.
"""
NumericPathProperties(coeff::Number) = NumericPathProperties(float(coeff), 0, 0, 0)

"""
Pretty print for NumericPathProperties
"""
Base.show(io::IO, pth::NumericPathProperties) = print(io, "NumericPathProperties($(pth.coeff), nsins=$(pth.nsins), ncos=$(pth.ncos), freq=$(pth.freq))")

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
    return PauliString(pstr.nqubits, pstr.operator, PathPropertiesType(pstr.coeff))
end

"""
    wrapcoefficients(psum::PauliSum)

Wrap the coefficients of a `PauliSum` into `NumericPathProperties`.
"""
function wrapcoefficients(psum::PauliSum)
    # by default wrap into `NumericPathProperties`
    return wrapcoefficients(psum, NumericPathProperties)
end

"""
    wrapcoefficients(psum::PauliSum, PathPropertiesType::Type{PP}) where {PP<:PathProperties}

Wrap the coefficients of a `PauliSum` into a custom `PathProperties` type.
For anything that is not natively supported by the library, you can subtype `PathProperties`.
A one-argument constructor of the custom `PathProperties` type from a coefficient must be defined.
"""
function wrapcoefficients(psum::PauliSum, PathPropertiesType::Type{PP}) where {PP<:PathProperties}
    return PauliSum(psum.nqubits, Dict(op => PathPropertiesType(coeff) for (op, coeff) in psum.op_dict))
end