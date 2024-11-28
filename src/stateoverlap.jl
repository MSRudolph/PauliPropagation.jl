### This file contains the functions to calculate the overlap between between backpropagated operators and the initial state.

"""
    overlapbyorthogonality(psum::PauliSum, orthogonalfunc::Function)

Overlap a `PauliSum` with a state or operator via function that returns true if a Pauli string is orthogonal and hence doesn't contribute.
An example `orthogonalfunc` is `containsXorY` which returns true if a Pauli string contains an X or Y Pauli.
If not orthogonal, then a Pauli string contributes with its coefficient.
This is particularly useful for overlaps with stabilizer states.
"""
function overlapbyorthogonality(psum::PauliSum, orthogonalfunc::Function)
    return overlapbyorthogonality(psum.op_dict, orthogonalfunc)
end

"""
    overlapbyorthogonality(pstr::PauliString, orthogonalfunc::Function)

Overlap a `PauliString` with a state or operator via function that returns true if the `PauliString` is orthogonal and hence has overlap 0.
 An example `orthogonalfunc` is `containsXorY` which returns true if the `PauliString` contains an X or Y Pauli.
If not orthogonal, then the overlap is the coefficient of the `PauliString`.
This is particularly useful for overlaps with stabilizer states.
"""
function overlapbyorthogonality(pstr::PauliString, orthogonalfunc::Function)
    return !orthogonalfunc(operator) * getnumcoeff(pstr.coeff)
end

"""
    overlapbyorthogonality(psum::Dict, orthogonalfunc::Function)

Overlap a Pauli sum dict with a state or operator via function that returns true if a Pauli string is orthogonal and hence doesn't contribute.
An example `orthogonalfunc` is `containsXorY` which returns true if a Pauli string contains an X or Y Pauli.
If not orthogonal, then a Pauli string contributes with its coefficient.
This is particularly useful for overlaps with stabilizer states.
"""
function overlapbyorthogonality(psum::Dict, orthogonalfunc::Function)
    val = keytype(psum)(0)
    for (operator, coeff) in psum
        if overlapbyorthogonality(operator, orthogonalfunc)
            val += getnumcoeff(coeff)
        end
    end
    return val
end

"""
    overlapbyorthogonality(pstr::PauliString, orthogonalfunc::Function)

Overlap an integer Pauli string with a state or operator via function that returns true if the Pauli string is orthogonal and hence has overlap 0.
 An example `orthogonalfunc` is `containsXorY` which returns true if the Pauli string contains an X or Y Pauli.
If not orthogonal, then the overlap is 1.
This is particularly useful for overlaps with stabilizer states.
"""
function overlapbyorthogonality(pstr::PauliStringType, orthogonalfunc::Function)
    return !orthogonalfunc(pstr)
end

## For the typical |0> or |+> cases
"""
    overlapwithzero(psum) 

Calculates the overlap of a Pauli sum with the zero state |0><0|
"""
overlapwithzero(psum) = overlapbyorthogonality(psum, orthogonaltozero)

"""
    overlapwithzero(psum) 

Calculates the overlap of a Pauli string with the zero state |0><0|
"""
orthogonaltozero(pstr) = containsXorY(pstr)

"""
    overlapwithplus(psum) 

Calculates the overlap of a Pauli sum with the plus state |+><+|
"""
overlapwithplus(psum) = overlapbyorthogonality(psum, orthogonaltoplus)

"""
    orthogonaltoplus(psum) 

Calculates the overlap of a Pauli string with the plus state |+><+|
"""
orthogonaltoplus(pstr) = containsYorZ(pstr)

# eval against |±i> not implemented


"""
    overlapwithmaxmixed(psum::PauliSum)

Calculates the overlap of a `PauliSum` with the maximally mixed state 1/2^n I.
"""
function overlapwithmaxmixed(psum::PauliSum)
    return overlapwithmaxmixed(psum.op_dict)
end

"""
    overlapwithmaxmixed(psum::Dict)

Calculates the overlap of a Pauli sum dict with the maximally mixed state 1/2^n I.
"""
function overlapwithmaxmixed(psum::Dict)
    IntType = keytype(psum)
    return get(psum, IntType(0), 0.0)
end

"""
    overlapwithpaulisum(psum1::PauliSum, psum2::PauliSum)

Calculates the overlap between two `PauliSum`s.
"""
function overlapwithpaulisum(psum1::PauliSum, psum2::PauliSum)
    return overlapwithpaulisum(psum1.op_dict, psum2.op_dict)
end

"""
    overlapwithpaulisum(psum1::PauliSum, psum2::PauliSum)

Calculates the overlap between two Pauli sum dicts.
"""
function overlapwithpaulisum(psum1::Dict, psum2::Dict)
    val = 0.0

    longer_psum = psum1
    shorter_psum = psum2

    # swap dicts around if op_dict is sparser
    if length(longer_psum) < length(shorter_psum)
        longer_psum, shorter_psum = shorter_psum, longer_psum
    end

    # looping over d2 (default initstate_dict) because we know that this one is sparser
    for operator in keys(shorter_psum)
        val += getnumcoeff(get(longer_psum, operator, 0.0)) * getnumcoeff(get(shorter_psum, operator, 0.0))
    end
    return val
end


"""
    filter(psum::PauliSum, filterfunc::Function)

Return a filtered `PauliSum` by removing all Pauli strings that satisfy the `filterfunc`.
"""
function filter(psum::PauliSum, filterfunc::Function)
    op_dict = filter(psum.op_dict, filterfunc)
    return PauliSum(psum.nqubits, op_dict)
end

"""
    filter!(psum::PauliSum, filterfunc::Function)

Filter a `PauliSum` in-place by removing all Pauli strings that satisfy the `filterfunc`.
"""
function filter!(psum::PauliSum, filterfunc::Function)
    filter!(psum.op_dict, filterfunc)
    return psum
end

"""
    filter(psum::Dict, filterfunc::Function)

Return a filtered Pauli sum dict by removing all Pauli strings that satisfy the `filterfunc`.
"""
function filter(psum::Dict, filterfunc::Function)
    return Dict(k => v for (k, v) in psum if !filterfunc(k))
end

"""
    filter!(psum::Dict, filterfunc::Function)

Filter a Pauli sum dict in-place by removing all Pauli strings that satisfy the `filterfunc`.
"""
function filter!(psum::Dict, filterfunc::Function)
    for pstr in keys(psum)
        if filterfunc(pstr)
            delete!(psum, pstr)
        end
    end
    return psum
end


# returns a new filtered dictionary, but doesn't overlap with anything
"""
    zerofilter(psum)

Return a filtered Pauli sum with only Pauli strings that are not orthogonal to the zero state |0><0|.
"""
zerofilter(psum) = filter(psum, containsXorY)

"""
    zerofilter!(psum)

Filter a Pauli sum in-place with only Pauli strings that are not orthogonal to the zero state |0><0|.
"""
zerofilter!(psum) = filter!(psum, containsXorY)

"""
    plusfilter(psum)

Return a filtered Pauli sum with only Pauli strings that are not orthogonal to the plus state |+><+|.
"""
plusfilter(psum) = filter(psum, containsYorZ)

"""
    zerofilter!(psum)

Filter a Pauli sum in-place with only Pauli strings that are not orthogonal to the plus state |+><+|.
"""
plusfilter!(psum) = filter!(psum, containsYorZ)


## Interface functions for extracting the numerical coefficients
"""
    getnumcoeff(val::Number)

Trivial function returning a numerical value of a number.
"""
function getnumcoeff(val::Number)
    return val
end

"""
    getnumcoeff(val::PathProperties)

Get the numerical coefficient of a `PathProperties` wrapper.
"""
function getnumcoeff(val::PathProperties)
    return val.coeff
end




