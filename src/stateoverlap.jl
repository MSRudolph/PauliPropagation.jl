### This file contains the functions to calculate the overlap between between backpropagated operators and the initial state.

## Evaluate with a rule 
function evalwithfilterfunction(op_dict, filterfunc)
    val = 0.0
    for (symbs, coeff) in op_dict
        contributes = !filterfunc(symbs)
        if contributes
            val += getnumcoeff(coeff)
        end
    end
    return val
end

## For the typical |0> case
evalagainstzero(op_dict) = evalwithfilterfunction(op_dict, containsXorY)
annihilatesatzero(op) = containsXorY(op)

evalagainstplus(op_dict) = evalwithfilterfunction(op_dict, containsYorZ)
annihilatesatplus(op) = containsYorZ(op)

# eval against |±i> not implemented

## Filter backpropagated operators
function filteroperators(op_dict, filterfunc)
    return Dict(k => v for (k, v) in op_dict if !filterfunc(k))
end

zerofilter(op_dict) = filteroperators(op_dict, containsXorY)

## Evaluate against initial state in dict form

function evalagainstdict(op_dict, initstate_dict)
    val = 0.0

    d1 = op_dict
    d2 = initstate_dict

    # swap dicts around if op_dict is sparser
    if length(d1) < length(d2)
        d1, d2 = d2, d1
    end

    # looping over d2 (default initstate_dict) because we know that this one is sparser
    for operator in keys(d2)
        val += getnum(get(d1, operator, 0.0)) * getnum(get(d2, operator, 0.0))
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




