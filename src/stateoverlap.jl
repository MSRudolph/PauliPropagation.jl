### This file contains the functions to calculate the overlap between between backpropagated pauli strings and states or general operators.

"""
    overlapbyorthogonality(orthogonalfunc::Function, psum::PauliSum)
    overlapbyorthogonality(orthogonalfunc::Function, pstr::PauliString)
    overlapbyorthogonality(orthogonalfunc::Function, pstr::Integer)

Overlap a `PauliSum`, `PauliString` or integer Pauli string with a state or operator via function 
that returns true if a Pauli string is orthogonal or not.
If a Pauli string is orthogonal, it does not contribute, otherwise it contributes with its coefficient.
This is particularly useful for overlaps with stabilizer states.
An example `orthogonalfunc` is `containsXorY` which returns true if a Pauli string contains an X or Y Pauli.
"""
function overlapbyorthogonality(orthogonalfunc::F, psum::PauliSum) where {F<:Function}
    if length(psum) == 0
        return 0.0
    end

    val = zero(numcoefftype(psum))
    for (pstr, coeff) in psum
        if overlapbyorthogonality(orthogonalfunc, pstr)
            val += tonumber(coeff)
        end
    end
    return val
end


function overlapbyorthogonality(orthogonalfunc::F, pstr::PauliString) where {F<:Function}
    return !orthogonalfunc(pstr) * tonumber(pstr.coeff)
end


function overlapbyorthogonality(orthogonalfunc::F, pstr::PauliStringType) where {F<:Function}
    return !orthogonalfunc(pstr)
end

## For the typical |0> or |+> cases
"""
    overlapwithzero(psum) 
    overlapwithzero(pstr)

Calculates the overlap of a Pauli sum with the zero state |0><0|,
i.e., Tr[psum * |0><0|] = <0|psum|0> or Tr[pstr * |0><0|] = <0|pstr|0>.
"""
overlapwithzero(pobj) = overlapbyorthogonality(containsXorY, pobj)

"""
    overlapwithplus(psum) 
    overlapwithplus(pstr)

Calculates the overlap of a Pauli sum or Pauli string with the plus state |+><+|,
i.e. Tr[psum * |+><+|] = <+|psum|+> or Tr[pstr * |+><+|] = <+|pstr|+>.
"""
overlapwithplus(psum) = overlapbyorthogonality(containsYorZ, psum)


# eval against |±i> not implemented


"""
    overlapwithcomputational(psum::PauliSum, onebitinds)
    overlapwithcomputational(pstr::PauliString, onebitinds)

Calculates the overlap of a Pauli sum or Pauli string with the computational basis state 
which has one-bits at all specified `indices` and zero-bits elsewhere.
If |x><x| is a computational basis state, it we compute Tr[psum * |x><x|] = <x|psum|x> or Tr[pstr * |x><x|] = <x|pstr|x>.
For example, `overlapwithcomputational(psum, [1,2,4])` returns the overlap with `|1101000...>`.
"""
function overlapwithcomputational(psum::PauliSum, onebitinds)
    val = zero(numcoefftype(psum))
    for (pstr, coeff) in psum
        val += tonumber(coeff) * _calcsignwithones(pstr, onebitinds)
    end
    return val
end


function overlapwithcomputational(pstr::PauliString, onebitinds)
    return _calcsignwithones(pstr.term, onebitinds) * tonumber(pstr.coeff)
end


function _calcsignwithones(pstr::PauliStringType, onebitinds)

    # factor is zero unless pstr is entirely I and Z
    if containsXorY(pstr)
        return 0
    end

    # factor is +-1 per the parity of pstr's Z=3 values at the bit=1 indices
    return (-1)^count(i -> getpauli(pstr, i) == 3, onebitinds)
end

"""
    overlapwithmaxmixed(psum::PauliSum)

Calculates the overlap of a `PauliSum` with the maximally mixed state I/2^n,
i.e., Tr[psum * I/2^n].
"""
function overlapwithmaxmixed(psum::PauliSum{TT,CT}) where {TT,CT}
    if length(psum) == 0
        return 0.0
    end

    NumType = numcoefftype(psum)

    return get(psum.terms, identitypauli(TT), zero(NumType))
end


"""
    scalarproduct(psum1::PauliSum, psum2::PauliSum)
    scalarproduct(pstr::PauliString, psum::PauliSum)
    scalarproduct(psum::PauliSum, pstr::PauliString)
    scalarproduct(pstr1::PauliString, pstr2::PauliString)

Calculates the scalar product between any combination of `PauliSum` and `PauliString`.
This  calculates the sum of the products of their coefficients for all Pauli strings that are present .
Important: This is not equivalent to the trace `Tr[psum1 * psum2]` but instead  `Tr[psum1 * psum2]/2^n`,
and equivalently for Pauli strings.
"""
function scalarproduct(psum1::PauliSum, psum2::PauliSum)
    if length(psum1) == 0 || length(psum2) == 0
        return 0.0
    end

    NumberType = numcoefftype(psum1)

    val = float(zero(NumberType))

    longer_psum = psum1.terms
    shorter_psum = psum2.terms

    # swap psums around if the other one is sparser
    if length(longer_psum) < length(shorter_psum)
        longer_psum, shorter_psum = shorter_psum, longer_psum
    end

    # looping over the shorter psum because we are only looking for collisions
    for pstr in keys(shorter_psum)
        val += tonumber(get(longer_psum, pstr, zero(NumberType))) * tonumber(get(shorter_psum, pstr, zero(NumberType)))
    end
    return val

end


function scalarproduct(pstr::PauliString, psum::PauliSum)
    return tonumber(getcoeff(psum, pstr)) * tonumber(pstr.coeff)

end

scalarproduct(psum::PauliSum, pstr::PauliString) = scalarproduct(pstr, psum)



function scalarproduct(pstr1::PauliString, pstr2::PauliString)
    if pstr1.term != pstr2.term
        return zero(numcoefftype(pstr1))
    end

    # at this point the two Pauli strings are the same
    # so we can just multiply the coefficients
    return tonumber(pstr1.coeff) * tonumber(pstr2.coeff)
end


"""
    filter(filterfunc::Function, psum::PauliSum)

Return a filtered `PauliSum` by removing all Pauli strings for which `filterfunc(pstr, coeff)` returns `false`.
"""
function Base.filter(filterfunc::F, psum::PauliSum) where {F<:Function}
    # iterating over dictionaries returns pairs like key=>value
    # so we need to unpack them to use the in-built Julia filter function
    filtered_terms = Base.filter(pair -> filterfunc(pair...), psum.terms)
    return PauliSum(psum.nqubits, filtered_terms)
end

"""
    filter!(filterfunc::Function, psum::PauliSum)

Filter a `PauliSum` in-place by removing all Pauli strings for which `filterfunc(pstr, coeff)` returns `false`.
"""
function Base.filter!(filterfunc::F, psum::PauliSum) where {F<:Function}
    # iterating over dictionaries returns pairs like key=>value
    # so we need to unpack them to use the in-built Julia filter function
    Base.filter!(pair -> filterfunc(pair...), psum.terms)
    return psum
end


# returns a new filtered dictionary, but doesn't overlap with anything
"""
    zerofilter(psum)

Return a filtered Pauli sum with only Pauli strings that are not orthogonal to the zero state |0><0|.
"""
zerofilter(psum) = filter((pstr, coeff) -> !containsXorY(pstr), psum)

"""
    zerofilter!(psum)

Filter a Pauli sum in-place with only Pauli strings that are not orthogonal to the zero state |0><0|.
"""
zerofilter!(psum) = filter!((pstr, coeff) -> !containsXorY(pstr), psum)

"""
    plusfilter(psum)

Return a filtered Pauli sum with only Pauli strings that are not orthogonal to the plus state |+><+|.
"""
plusfilter(psum) = filter((pstr, coeff) -> !containsYorZ(pstr), psum)

"""
    zerofilter!(psum)

Filter a Pauli sum in-place with only Pauli strings that are not orthogonal to the plus state |+><+|.
"""
plusfilter!(psum) = filter!((pstr, coeff) -> !containsYorZ(pstr), psum)




