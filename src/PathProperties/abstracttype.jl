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
Pretty print for PathProperties
"""
function Base.show(io::IO, pth::PProp) where {PProp<:PathProperties}
    print(io, "$PProp(")
    for (i, field) in enumerate(fieldnames(PProp))
        print(io, "$(field)=$(getfield(pth, field))")
        if i < length(fieldnames(PProp))
            print(io, ", ")
        end
    end
    print(io, ")")

end

import Base: *
"""
Multiplication of the `coeff` field in a `PathProperties` object with a number.
Requires that the `PathProperties` object has a `coeff` field defined which will be multiplied.
"""
function *(path::PProp, val::Number) where {PProp<:PathProperties}
    # multiply the coefficient on the `coeff` field with the value and leave the rest unchanged.
    fields = fieldnames(PProp)

    if :coeff ∉ fields
        throw("The $(PProp) object does not have a field `coeff` to use the `*` operation.")
    end

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
Requires that the `PathProperties` object has a `coeff` field defined which will be multiplied.
"""
function *(val::Number, path::PProp) where {PProp<:PathProperties}
    return path * val
end

import Base: +
"""
Addition of two `PathProperties` objects of equal concrete type.
Adds the `coeff` fields and takes the minimum of the other fields.
Requires that the `PathProperties` object has a `coeff` field defined.
"""
function +(path1::PProp, path2::PProp) where {PProp<:PathProperties}
    fields = fieldnames(PProp)

    if :coeff ∉ fields
        throw("The $(PProp) object does not have a field `coeff` to use the `+` operation.")
    end

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

import Base: float
"""
    float(path::PathProperties)

Returns an equivalent `PathProperties` object where float() is applied to the `coeff` field.
"""
function float(path::PProp) where {PProp<:PathProperties}
    fields = fieldnames(PProp)

    if :coeff ∉ fields
        throw("The $(PProp) object does not have a field `coeff` to use the `float` operation.")
    end

    # update the `coeff` only
    function updateval(fval, fname)
        if fname == :coeff
            fval = float(fval)
        end
        return fval
    end

    return PProp((updateval(getfield(path, fname), fname) for fname in fields)...)
end


"""
    tonumber(path::PathProperties)

Get the numerical coefficient of a `PathProperties` wrapper.
"""
function tonumber(path::PProp) where {PProp<:PathProperties}
    if !hasfield(PProp, :coeff)
        throw("The $(PProp) object does not have a field `coeff`.
        Consider defining a `tonumber(path::$(PProp))` method.")
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