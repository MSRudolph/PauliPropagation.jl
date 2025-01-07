### abstracttype.jl
##
# PathProperties abstract type and methods.
# Subtype PathProperties to wrap coefficients and record custom path properties.
# It is important for subtyping structs to have the `coeff` field where current numerical coefficient is stored.
# Additionally, it is assumed that, when paths are merged, the `coeff` fields are added 
# and the other fields are taken as the minimum between the respective fields on the two paths.
# This behavior can be changed but then more methods need to be defined.
##
###


"""
Abstract type for wrapping coefficients and record custom path properties
"""
abstract type PathProperties end

"""
PathProperties are expected to be immutable.
"""
Base.copy(path::PathProperties) = path

"""
Pretty print for PathProperties
"""
function Base.show(io::IO, pth::PProp) where {PProp<:PathProperties}
    print_string = "$PProp("
    for (i, field) in enumerate(fieldnames(PProp))
        print_string *= "$(field)=$(getfield(pth, field))"
        if i < length(fieldnames(PProp))
            print_string *= ", "
        end
    end
    print_string *= ")"

    print(io, print_string)

end

import Base: *
"""
Multiplication of the `coeff` field in a `PathProperties` object with a number.
Requires that the `PathProperties` object has a `coeff` field as the first field.
"""
function *(path::PProp, val::Number) where {PProp<:PathProperties}
    # multiply the coefficient on the `coeff` field with the value and leave the rest unchanged.
    fields = fieldnames(PProp)

    # update the `coeff` only
    function updateval(fval, fname)
        if fname == :coeff
            fval *= val
        end
        return fval
    end

    return PProp((updateval(getfield(path, fname), fname) for fname in fields)...)
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

    # add `coeff` fields and take the minimum of the other fields
    function updateval(fval1, fval2, fname)
        if fname == :coeff
            return fval1 + fval2
        else
            return min(fval1, fval2)
        end
    end

    return PProp((updateval(getfield(path1, fname), getfield(path2, fname), fname) for fname in fields)...)
end

"""
    numcoefftype(path::PathProperties)

Return the type of the coefficient `coeff` in a `PathProperties` object if applicable.
"""
function numcoefftype(path::PProp) where {PProp<:PathProperties}
    return coefftype(path)
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
    wrapcoefficients(pstr::PauliString, PathPropertiesType<:PathProperties)

Wrap the coefficient of a `PauliString` into a custom `PathProperties` type.
For anything that is not natively supported by the library, you can subtype `PathProperties`.
A one-argument constructor of the custom `PathProperties` type from a coefficient must be defined.
"""
function wrapcoefficients(pstr::PauliString, ::Type{PProp}) where {PProp<:PathProperties}
    # the one-argument constructor of your PathProperties type must be defined
    pprop = try
        PProp(pstr.coeff)
    catch MethodError
        throw(
            "The constructor `$(PProp)(coeff)` is not defined for the $(PProp) type. " *
            "Either construct a PauliString with wrapped coefficient or define the `$(PProp)(coeff)` constructor."
        )
    end

    return PauliString(pstr.nqubits, pstr.term, pprop)
end

"""
    wrapcoefficients(psum::PauliSum, PathPropertiesType<:PathProperties)

Wrap the coefficients of a `PauliSum` into a custom `PathProperties` type.
For anything that is not natively supported by the library, you can subtype `PathProperties`.
A one-argument constructor of the custom `PathProperties` type from a coefficient must be defined.
"""
function wrapcoefficients(psum::PauliSum, ::Type{PProp}) where {PProp<:PathProperties}
    if length(psum) == 0
        throw("The PauliSum is empty.")
    end

    try
        _, dummy_coeff = first(psum.terms)
        PProp(dummy_coeff)
    catch MethodError
        throw(
            "The constructor `$(PProp)(coeff)` is not defined for the $(PProp) type. " *
            "Either construct a PauliSum with wrapped coefficient or define the `$(PProp)(coeff)` constructor.")
    end

    return PauliSum(psum.nqubits, Dict(pstr => PProp(coeff) for (pstr, coeff) in psum.terms))
end