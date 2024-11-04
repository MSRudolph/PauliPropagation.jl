using Statistics



"""
Function to estimate the average error of a truncated circuit simulation using Monte Carlo sampling.
Returns the mean squared error of the truncated Pauli propagation simulation averaged over all parameters.

The function takes a circuit, a PauliString, and the number of Monte Carlo samples. Thetas will be set to ones and split_probabilities to 0.5 by default.
Importantly, the `kwargs` can be used to set the truncation parameters of the Pauli propagation.
"""
function estimateaverageerror(circ, pstr::PauliString, n_mcsamples::Integer; circuit_reversed=false, kwargs...)
    # length(thetas) should be equal to the number of parametrized gates in the circuit
    thetas = ones(countparameters(circ))
    # length(split_probabilities) should be equal to the number of total gates in the circuit
    split_probabilities = ones(length(circ)) .* 0.5

    outcome_array = zeros(n_mcsamples)

    # TODO: Need to provide initial state overlapfunction

    return estimateaverageerror!(circ, pstr, outcome_array, thetas, split_probabilities; circuit_reversed=circuit_reversed, kwargs...)

end



"""
Function to estimate the average error of a truncated circuit simulation using Monte Carlo sampling.
Returns the mean squared error of the truncated Pauli propagation simulation averaged over `thetas`.

The function takes a circuit, a PauliString, the number of Monte Carlo samples, an array of thetas, and an array of split_probabilities. 
Importantly, the `kwargs` can be used to set the truncation parameters of the Pauli propagation.

The length the `thetas` vector should be equal to the number of parametrized gates in the circuit. This can be awkward, for example for Pauli gates, 
because they have parameters but in the average error estimation they are not used. But the thetas matter for non-splitting parametrized gates.

The the length of the `split_probabilities` vector should be equal to the number of total gates in the circuit. 
This can be awkward, for example for non-splitting gates. But these matter for the splitting gates.
"""
function estimateaverageerror(circ, pstr::PauliString, n_mcsamples::Integer, thetas::AbstractArray, split_probabilities::AbstractArray; circuit_reversed=false, kwargs...)
    # length(thetas) should be equal to the number of parametrized gates in the circuit
    @assert length(thetas) == countparameters(circ)
    # length(split_probabilities) should be equal to the number of total gates in the circuit
    @assert length(split_probabilities) == length(circ)

    outcome_array = zeros(n_mcsamples)

    # TODO: Need to provide initial state overlapfunction

    return estimateaverageerror!(circ, pstr, outcome_array, thetas, split_probabilities; circuit_reversed=circuit_reversed, kwargs...)

end


"""
In-place version of `estimateaverageerror`. This function takes an array of length m_mcsamples as an argument and modifies it in-place.
"""
function estimateaverageerror!(circ, pstr::PauliString, outcome_array::AbstractVector, thetas::AbstractArray, split_probabilities::AbstractArray; circuit_reversed=false, kwargs...)
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
        # if not valid, then the error is coeff * overlapwithzero(final_operator), if valid, then the error is 0
        outcome_array[ii] = final_pstr.coeff * overlapwithzero(final_pstr.operator) * (1 - is_valid)
    end

    return mean(outcome_array)

end

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

        coeff, pauli = mcapply(gate, pauli, coeff, thetas[param_idx], split_probabilities[prob_idx]; kwargs...)

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
    return PauliString(pstr.nqubits, pauli, coeff), is_valid
end

# if we don't define an mcapply() function, assume that it doesn't split and that the apply function does the job
mcapply(gate, pauli, coeff, theta, args...) = apply(gate, pauli, theta, coeff) # TODO: have more sensible default functions once the numerical certificate is figured out.

# mcapply for Pauli gates
function mcapply(gate::PauliGateUnion, pauli, coeff, theta, split_prob=0.5)  # 

    if !commutes(gate, pauli)
        # if the gate does not commute with the operator, split with probability split_prob int a random direction.
        if rand() > split_prob
            # TODO: make sure this works in general and not just for numerical coefficients. For example, add to freq and nsins.
            sign, pauli = getnewoperator(gate, pauli)
        end
    end
    # the sign has abs(sign) = 1, so updating the coefficient with abs(sign) doesn't do anything
    return coeff, pauli
end

function mcapply(gate::AmplitudeDampingNoise, args...)  # 
    # TODO
    println("AmplitudeDampingNoise not implemented yet")
end
