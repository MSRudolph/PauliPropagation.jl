### mcgenerics.jl
##
# This file contains the foundational functions for Monte Carlo sampling with respect to 2-norms of coefficients. 
##
###


# `mcapply()` function for a `CliffordGate` is just the `apply()` function because it does not split.
mcapply(gate::CliffordGate, pstr::PauliString, theta, split_probability, rng) = PauliString(pstr.nqubits, apply(gate, pstr.term, pstr.coeff)...)


# MC apply function for a `MaskedPauliRotation`.
# This will error if a `PauliRotation` is not converted to a `MaskedPauliRotation` before calling this function.
# If the gate commutes with the pauli string, the pauli string is left unchanged. 
# Else the pauli string is split off with a probability 1 - `split_prob`.
function mcapply(gate::MaskedPauliRotation, pstr::PauliString, theta, split_prob, rng)

    if commutes(gate, pstr.term)
        return pstr
    end
    coeff = pstr.coeff

    # if the gate does not commute with the pauli string, remain with probability `split_prob` and split off with probability 1 - `split_prob`.
    if rand(rng) < split_prob
        new_term = pstr.term
        coeff *= cos(theta)
        # Pauli doesn't get changed      
    else
        # branch into the new Pauli string, ifnore ignore the sign
        new_term, sign = getnewpaulistring(gate, pstr.term)
        coeff *= sin(theta) * sign
    end

    # the sign has abs(sign) = 1, so updating the coefficient with abs(sign) doesn't do anything
    return PauliString(pstr.nqubits, new_term, coeff)
end

# Perform a single Monte Carlo propagation of a Pauli string through an already reversed circuit. Returns the final Pauli string and a boolean indicating whether the path was truncated.
# It further assumes that the `thetas` and `split_probabilities` are already correctly calculated and provided as arguments. 
# In general they will be vectors, but they can also be real numbers.
function montecarlopropagation(circ, pstr::PauliString, thetas, split_probabilities, rng; min_coeff=0., max_weight=Inf, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing)

    param_idx = length(thetas)
    prob_idx = length(split_probabilities)

    for gate in circ
        # apply the gate to the Pauli string
        # if the gate splits, the Pauli string is split with a probability 1 - split_prob
        pstr = mcapply(gate, pstr, _getelmt(thetas, param_idx), _getelmt(split_probabilities, prob_idx), rng)

        # decrement the parameter index if gate is parametrized
        if isa(gate, ParametrizedGate) && param_idx > 0
            param_idx -= 1
        end
        # always decrement the probability index
        # will only be zero at the end of the circuit
        prob_idx -= 1
    end
    return pstr
end

function _calculatesplitprobabilities(circ::AbstractArray, thetas::AbstractArray)
    if length(thetas) != countparameters(circ)
        throw("Vector `thetas` must have same length the number of parametrized gates in `circ`.")
    end

    split_probabilities = zeros(length(circ))

    theta_idx = 1
    for (ii, gate) in enumerate(circ)
        if isa(gate, ParametrizedGate)
            split_probabilities[ii] = _calculatesplitprobabilities(gate, thetas[theta_idx])
            theta_idx += 1
        end
    end
    return split_probabilities
end


# Function that automatically calculates the splitting probability of the gates in the circuit based on a one number theta.
# This assumes that the circuit consists only of `PauliRotation` -`CliffordGate`.

_calculatesplitprobabilities(gate::PauliRotationUnion, theta::Number) = cos(theta)^2

# helper function to get element from array or return the number itself
_getelmt(arr::AbstractArray, idx::Integer) = arr[idx]
_getelmt(num::Number, idx::Integer) = num

