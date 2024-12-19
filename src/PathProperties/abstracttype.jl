### abstracttype.jl
##
# PathProperties abstract type and methods.
# Subtype PathProperties to wrap coefficients and record custom path properties.
# Properties to track could for example be the number of times a path branched at a PauliRotation gate.
# Another could be the so-called path weight, which is the sum of Pauli weights of the operator after each gate.
# It is important for subtyping structs to have the `coeff` field as the first field where current numerical coefficient is stored.
# If the struct is parametrized, the first parameter should be the type of the coefficient.
# This behavior can be changed but then more methods need to be defined.
##
###


"""
Abstract type for wrapping coefficients and record custom path properties
"""
abstract type PathProperties end

"""
Pretty print for PathProperties
"""
Base.show(io::IO, path::PathProperties) = print(io, "PathProperties($(typeof(path.coeff)), nsins=$(path.nsins), ncos=$(path.ncos), freq=$(path.freq))")

"""
PathProperties are expected to be immutable.
"""
Base.copy(path::PathProperties) = path

import Base: *
"""
Multiplication of the `coeff` field in a `PathProperties` object with a number.
Requires that the `PathProperties` object has a `coeff` field as the first field.
"""
function *(path::PProp, val::Number) where {PProp<:PathProperties}
    # multiply the coefficient on the `coeff` field with the value and leave the rest unchanged.
    fields = fieldnames(PProp)
    return PProp(getfield(path, :coeff) * val, (getfield(path, fields[ii]) for ii in 2:length(fields))...)
end

"""
Multiplication of a `PathProperties` object with a number.
Requires that the `PathProperties` object has a `coeff` field as the first field which will be multiplied.
"""
function *(val::Number, path::PProp) where {PProp<:PathProperties}
    return path * val
end

import Base: +
"""
Addition of two `PathProperties` objects of equal concrete type.
Adds the `coeff` fields and takes the minimum of the other fields.
Requires that the `PathProperties` object has a `coeff` field as the first field.
"""
function +(path1::PProp, path2::PProp) where {PProp<:PathProperties}
    fields = fieldnames(PProp)
    return PProp(getfield(path1, :coeff) + getfield(path2, :coeff), (min(getfield(path1, fields[ii]), getfield(path2, fields[ii])) for ii in 2:length(fields))...)
end

"""
    numcoefftype(path::PathProperties)

Return the type of the coefficient `coeff` in a `PathProperties` object if applicable.
"""
function numcoefftype(path::PProp) where {PProp<:PathProperties}
    return coefftype(path)
end

"""
    numcoefftype(::Type{PathProperties})

Return the first parameter in a in a parametrized `PathProperties` object if applicable.
"""
function numcoefftype(::Type{PProp}) where {PProp<:PathProperties}
    return coefftype(PProp)
end

"""
    coefftype(path::PathProperties)

Return the type of the coefficient `coeff` in a `PathProperties` object if applicable.
"""
function coefftype(path::PProp) where {PProp<:PathProperties}
    if !hasfield(PProp, :coeff)
        throw("The $(PProp) object does not have a field `coeff` to determine the numerical coefficient type.
        Consider defining a `coefftype(path::$(PProp))` method.")
    end
    return typeof(path.coeff)
end

"""
    coefftype(::Type{PathProperties})

Return the first parameter in a in a parametrized `PathProperties` object if applicable.
"""
function coefftype(::Type{PProp}) where {PProp<:PathProperties}
    if length(PProp.parameters) == 0
        throw("The $(PProp) type is not parametrized to determine the numerical coefficient type." *
              "Consider defining a `coefftype(path::Type{$(PProp)})` method.")
    end

    CT = PProp.parameters[1]

    if CT <: Number
        return CT
    else
        throw("The first parameter of the $(PProp) type is $CT and not a number type." *
              "Consider defining a `coefftype(path::Type{$(PProp)})` method.")
    end

end

"""
    tonumber(val::PathProperties)

Get the numerical coefficient of a `PathProperties` wrapper.
"""
function tonumber(path::PProp) where {PProp<:PathProperties}
    if !hasfield(PProp, :coeff)
        throw("The $(PProp) object does not have a field `coeff` to determine the numerical coefficient type.
        Consider defining a `coefftype(path::$(PProp))` method.")
    end
    return path.coeff
end


"""
    wrapcoefficients(pstr::PauliString, PathPropertiesType::Type{PP}) where {PP<:PathProperties}

Wrap the coefficient of a `PauliString` into a custom `PathProperties` type.
For anything that is not natively supported by the library, you can subtype `PathProperties`.
A one-argument constructor of the custom `PathProperties` type from a coefficient must be defined.
"""
function wrapcoefficients(pstr::PauliString, ::Type{PProp}) where {PProp<:PathProperties}
    # the one-argument constructor of your PathProperties type must be defined
    return PauliString(pstr.nqubits, pstr.term, PProp(pstr.coeff))
end

"""
    wrapcoefficients(psum::PauliSum, PathPropertiesType::Type{PP}) where {PP<:PathProperties}

Wrap the coefficients of a `PauliSum` into a custom `PathProperties` type.
For anything that is not natively supported by the library, you can subtype `PathProperties`.
A one-argument constructor of the custom `PathProperties` type from a coefficient must be defined.
"""
function wrapcoefficients(psum::PauliSum, ::Type{PProp}) where {PProp<:PathProperties}
    return PauliSum(psum.nqubits, Dict(pstr => PProp(coeff) for (pstr, coeff) in psum.terms))
end