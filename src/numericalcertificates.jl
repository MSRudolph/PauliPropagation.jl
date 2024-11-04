using Statistics



"""
Function to estimate the average error of a truncated Pauli propagation simulation using Monte Carlo sampling.
This version of the function is only valid for a circuit consisting of `PauliGates` non-parametrized non-splitting gates (for example Cliffords).
Returns the mean squared error of the truncated Pauli propagation simulation averaged over the angles of all `PauliGate`s in the interval [-pi, pi].

An initial state overlap function `initialstatefunc` can be provided to calculate the overlap of the backpropagated Pauli operators with the initial state.
Importantly, the `kwargs` can be used to set the truncation parameters of the Pauli propagation.
"""
function estimateaverageerror(circ, pstr::PauliString, n_mcsamples::Integer; initialstatefunc=overlapwithzero, circuit_reversed=false, kwargs...)
    # this function is only valid for PauliGates and non-splitting non-parametrized gates
    # this will not not error for non-parametrized splitting gates, e.g. T-gates. 
    @assert all(g -> isa(g, PauliGateUnion) || isa(g, StaticGate), circ)

    # integration range over thetas is [-pi, pi]
    thetas = ones(countparameters(circ)) * pi

    return estimateaverageerror(circ, pstr, n_mcsamples, thetas; initialstatefunc=initialstatefunc, circuit_reversed=circuit_reversed, kwargs...)
end



"""
Function to estimate the average error of a truncated circuit simulation using Monte Carlo sampling.
Returns the mean squared error of the truncated Pauli propagation simulation averaged over the `thetas`∈ [theta, theta] of the angle `theta` of each `PauliGate`.

The length the `thetas` vector should be equal to the number of parametrized gates in the circuit. For `PauliGates`, the value should be the integration range of the parameters around zero.
For other currently supported parametrized gates, potential splitting probabilities can be derived from the parameters (e.g., for `AmplitudeDampingNoise`). 
We currently support no non-parametrized splitting gates. This may change in the future.

An initial state overlap function `initialstatefunc` can be provided to calculate the overlap of the backpropagated Pauli operators with the initial state.
Importantly, the `kwargs` can be used to set the truncation parameters of the Pauli propagation.
"""
function estimateaverageerror(circ, pstr::PauliString, n_mcsamples::Integer, thetas::AbstractArray; initialstatefunc=overlapwithzero, circuit_reversed=false, kwargs...)
    # length(thetas) should be equal to the number of parametrized gates in the circuit
    @assert length(thetas) == countparameters(circ)

    split_probabilities = _calculatesplitprobabilities(circ, thetas)

    outcome_array = zeros(n_mcsamples)

    return estimateaverageerror!(circ, pstr, outcome_array, thetas, split_probabilities; initialstatefunc=initialstatefunc, circuit_reversed=circuit_reversed, kwargs...)

end

"""
In-place version of `estimateaverageerror`. This function takes an array of length n_mcsamples as an argument and modifies it in-place. 
It further assumes that the `thetas` and `split_probabilities` are already correctly calculated and provided as arguments.
"""
function estimateaverageerror!(circ, pstr::PauliString, outcome_array::AbstractVector, thetas::AbstractArray, split_probabilities::AbstractArray; initialstatefunc=paulioverlapwithzero, circuit_reversed=false, kwargs...)
    # This function takes an outcome_array as an argument and modifies it in-place.

    # length(thetas) should be equal to the number of parametrized gates in the circuit
    @assert length(thetas) == count(g -> isa(g, ParametrizedGate), circ)
    # length(split_probabilities) should be equal to the number of total gates in the circuit
    @assert length(split_probabilities) == length(circ)

    # reverse the circuit once 
    if circuit_reversed
        circ = reverse(circ)
    end

    n_mcsamples = length(outcome_array)
    @threads for ii in 1:n_mcsamples
        final_pstr, is_valid = montecarlopropagation(circ, pstr, thetas, split_probabilities; circuit_reversed=true, kwargs...)

        # multiply the coefficient of the backpropagated Pauli with the overlap with the initial state (#TODO: provide initial state)
        # and then multply with (1 - is_valid) to get the final outcome.
        # if not valid, then the error is coeff * initial_state_func(final_operator), if valid, then the error is 0
        outcome_array[ii] = final_pstr.coeff * initialstatefunc(final_pstr.operator) * (1 - is_valid)
    end

    return mean(outcome_array)

end

## TODO: Have easier to use versions of the functions that don't require the user to calculate the split probabilities and thetas.
function montecarlopropagation(circ, pstr::PauliString, thetas::AbstractArray, split_probabilities::AbstractArray; circuit_reversed=false, max_weight=Inf, max_freq=Inf, max_sins=Inf, kwargs...)
    # Reverse the circ if it is not already done. Allocates memory.
    if circuit_reversed
        circ = reverse(circ)
    end
    param_idx = length(thetas)
    prob_idx = length(split_probabilities)

    is_valid = true

    pauli = pstr.operator
    coeff = pstr.coeff
    for gate in circ

        coeff, pauli = mcapply(gate, pauli, coeff, thetas[param_idx], split_probabilities[prob_idx]; kwargs...) # TODO: this currently allocates a lot of memory.

        # check truncations
        # TODO: make this properly
        if truncateweight(pauli, max_weight)
            is_valid = false
        elseif truncatefrequency(coeff, max_freq)
            is_valid = false
        elseif truncatesins(coeff, max_sins)
            is_valid = false
        end

        # TODO: Custom truncation functions?

        # decrement the parameter index if gate is parametrized
        if isa(gate, ParametrizedGate)
            param_idx -= 1
        end
        # always decrement the probability index
        prob_idx -= 1
    end
    return PauliString{typeof(pstr).parameters...}(pstr.nqubits, pauli, coeff), is_valid
end

# if we don't define an mcapply() function, assume that it doesn't split and that the apply function does the job
mcapply(gate, pauli, coeff, theta, split_probability; kwargs...) = apply(gate, pauli, theta, coeff; kwargs...) # TODO: have more sensible default functions once the numerical certificate is figured out.

## Monte Carlo apply functions
"""
MC apply function for Pauli gates. If the gate commutes with the operator, the operator is left unchanged. Else the operator is split off with a probability 1 - `split_prob`.
"""
function mcapply(gate::PauliGateUnion, pauli, coeff, theta, split_prob=0.5; kwargs...)  # 

    if !commutes(gate, pauli)
        # if the gate does not commute with the operator, remain with probability `split_prob` and split off with probability 1 - `split_prob`.
        if rand() > split_prob
            # TODO: make sure this works in general and not just for numerical coefficients. For example, add to freq and nsins.
            sign, pauli = getnewoperator(gate, pauli)
        end
    end
    # the sign has abs(sign) = 1, so updating the coefficient with abs(sign) doesn't do anything
    return coeff, pauli
end

function mcapply(gate::AmplitudeDampingNoise, args...; kwargs...)  # 
    # TODO
    throw("AmplitudeDampingNoise not implemented yet")
end



### Utilities
"""
Function that automatically calculates the splitting probabilities of the gates in the circuit based on a vector of thetas.
For Pauli gates, the theta value is interpreted as the limits of the integration [-theta, theta].
For AmplitudeDampingNoise, the splitting probability is the damping rate.
"""
function _calculatesplitprobabilities(circ, thetas)
    split_probabilities = zeros(length(circ))

    theta_idx = 1
    for (ii, gate) in enumerate(circ)
        if isa(gate, PauliGateUnion)
            # integration range
            r = thetas[theta_idx]
            split_probabilities[ii] = 0.5 * (1 + sin(2r) / (2r))
        elseif isa(gate, AmplitudeDampingNoise)
            # probability of splitting = damping rate. Need to rescale diagonal accordingly.
            split_probabilities[ii] = thetas[theta_idx]
        end
        if isa(gate, ParametrizedGate)
            theta_idx += 1
        end
    end
    return split_probabilities
end