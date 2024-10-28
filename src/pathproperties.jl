## This type can be used to wrap coefficients and record custom properties
abstract type PathProperties end

## Specific PathProperties
mutable struct NumericPathProperties <: PathProperties
    coeff::Float64
    nsins::Int
    ncos::Int
    freq::Int
end

## This 1-argument constructor needs to be defined for any PathProperties type
NumericPathProperties(coeff) = NumericPathProperties(coeff, 0, 0, 0)

import Base: *
function *(pth::PathProperties, val::Number)
    pth.coeff *= val
    return pth
end
import Base: copy
function copy(path_properties::PathProperties)
    return typeof(path_properties)(path_properties.coeff, path_properties.nsins, path_properties.ncos, path_properties.freq)
end

Base.show(io::IO, pth::PathProperties) = print(io, "PathProperties($(typeof(pth.coeff)), nsins=$(pth.nsins), ncos=$(pth.ncos), freq=$(pth.freq))")

Base.show(io::IO, pth::NumericPathProperties) = print(io, "NumericPathProperties($(pth.coeff), nsins=$(pth.nsins), ncos=$(pth.ncos), freq=$(pth.freq))")


## Wrapping PauliString and PauliSum in PathProperties
function wrapcoefficients(pstr::PauliString, PathPropertiesType::Type{PP}) where {PP<:PathProperties}
    # the one-argument constructor of your PathProperties type must be defined
    return PauliString(pstr.nqubits, pstr.operator, PathPropertiesType(pstr.coeff))
end

function wrapcoefficients(psum::PauliSum, PathPropertiesType::Type{PP}) where {PP<:PathProperties}
    return PauliSum(psum.nqubits, Dict(op => PathPropertiesType(coeff) for (op, coeff) in psum.op_dict))
end