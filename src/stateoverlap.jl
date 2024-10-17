### This file contains the functions to calculate the overlap between between backpropagated operators and the initial state.

## Evaluate with a rule 
function overlapbyorthogonality(psum::PauliSum, orthogonalfunc)
    return overlapbyorthogonality(psum.op_dict, orthogonalfunc)
end

function overlapbyorthogonality(pstr::PauliString, orthogonalfunc)
    return !orthogonalfunc(operator) * getnumcoeff(pstr.coeff)
end

function overlapbyorthogonality(op_dict::Dict, orthogonalfunc)
    val = 0.0
    for (operator, coeff) in op_dict
        if !orthogonalfunc(operator)
            val += getnumcoeff(coeff)
        end
    end
    return val
end

## For the typical |0> or |+> cases
overlapwithzero(op_dict) = overlapbyorthogonality(op_dict, orthogonaltozero)
orthogonaltozero(op) = containsXorY(op)

overlapwithplus(op_dict) = overlapbyorthogonality(op_dict, orthogonaltoplus)
orthogonaltoplus(op) = containsYorZ(op)

# eval against |±i> not implemented

# TODO: Implement with overlap maximally mixed once we have the operator interface

## Evaluate the overlap (or trace) between two PauliSum objects
function overlapwithpaulisum(psum1::PauliSum, psum2::PauliSum)
    return overlapwithdict(psum1.op_dict, psum2.op_dict)
end

function overlapwithdict(op_dict::Dict, initstate_dict::Dict)
    val = 0.0

    d1 = op_dict
    d2 = initstate_dict

    # swap dicts around if op_dict is sparser
    if length(d1) < length(d2)
        d1, d2 = d2, d1
    end

    # looping over d2 (default initstate_dict) because we know that this one is sparser
    for operator in keys(d2)
        val += getnumcoeff(get(d1, operator, 0.0)) * getnumcoeff(get(d2, operator, 0.0))
    end
    return val
end


## Filter backpropagated operators
function filterdict(op_dict::Dict, filterfunc)
    return Dict(k => v for (k, v) in op_dict if !filterfunc(k))
end

function filterpaulisum(psum::PauliSum, filterfunc)
    op_dict = filterdict(psum.op_dict, filterfunc)
    return PauliSum(psum.nqubits, op_dict)
end

# returns a new filtered dictionary, but doesn't overlap with anything
zerofilter(psum) = filterdict(psum, containsXorY)
plusfilter(psum) = filterdict(psum, containsYorZ)


## Interface functions for extracting the numerical coefficients
function getnumcoeff(val::Real)
    return val
end

function getnumcoeff(val::NumericPathProperties)
    return val.coeff
end




