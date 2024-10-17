### This file contains the functions to calculate the overlap between between backpropagated operators and the initial state.

## Evaluate with a rule 
function overlapbyorthogonality(op_dict, orthogonalfunc)
    val = 0.0
    for (operator, coeff) in op_dict
        if !orthogonalfunc(operator)
            val += getnumcoeff(coeff)
        end
    end
    return val
end

## For the typical |0> case
overlapwithzero(op_dict) = overlapbyorthogonality(op_dict, orthogonaltozero)
orthogonaltozero(op) = containsXorY(op)

overlapwithplus(op_dict) = overlapbyorthogonality(op_dict, orthogonaltoplus)
orthogonaltoplus(op) = containsYorZ(op)

# eval against |±i> not implemented

# TODO: Implement with overlap maximally mixed once we have the operator interface


## Filter backpropagated operators
function filterdict(op_dict, filterfunc)
    return Dict(k => v for (k, v) in op_dict if !filterfunc(k))
end

# returns a new filtered dictionary, but doesn't overlap with anything
zerofilter(op_dict) = filterdict(op_dict, containsXorY)
plusfilter(op_dict) = filterdict(op_dict, containsYorZ)

## Evaluate against initial state in dict form
function overlapwithdict(op_dict, initstate_dict)  # TODO: change name and functionality to operator interface
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


## Interface functions for extracting the numerical coefficients
function getnumcoeff(val::Real)
    return val
end

function getnumcoeff(val::NumericPathProperties)
    return val.coeff
end




