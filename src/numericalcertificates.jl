using Statistics

"""
Function to estimate the average error of a truncated circuit simulation using Monte Carlo sampling.
Returns the mean squared error of the truncated Pauli propagation simulation averaged over `thetas`.

The function takes a circuit, a PauliString, an array of thetas, an array of split_probabilities, and the number of Monte Carlo samples. 
Importantly, the `kwargs` can be used to set the truncation parameters of the Pauli propagation.

The length the `thetas` vector should be equal to the number of parametrized gates in the circuit. This can be awkward, for example for Pauli gates, 
because they have parameters but in the average error estimation they are not used. But the thetas matter for non-splitting parametrized gates.

The the length of the `split_probabilities` vector should be equal to the number of total gates in the circuit. 
This can be awkward, for example for non-splitting gates. But these matter for the splitting gates.
"""
function estimateaverageerror(circ, pstr::PauliString, thetas::AbstractArray, split_probabilities::AbstractArray, n_mcsamples::Integer; kwargs...)
    # length(thetas) should be equal to the number of parametrized gates in the circuit
    @assert length(thetas) == countparameters(circ)
    # length(split_probabilities) should be equal to the number of total gates in the circuit
    @assert length(split_probabilities) == length(circ)

    outcome_array = zeros(n_mcsamples)

    # TODO: Need to provide initial state overlapfunction

    return estimateaverageerror!(circ, pstr, thetas, split_probabilities, outcome_array; kwargs...)

end

function estimateaverageerror!(circ, pstr::PauliString, thetas::AbstractArray, split_probabilities::AbstractArray, outcome_array::AbstractVector; kwargs...)
    # This function takes an outcome_array as an argument and modifies it in-place.

    # length(thetas) should be equal to the number of parametrized gates in the circuit
    @assert length(thetas) == count(g -> isa(g, ParametrizedGate), circ)
    # length(split_probabilities) should be equal to the number of total gates in the circuit
    @assert length(split_probabilities) == length(circ)

    # reverse the circuit once 
    circ = reverse(circ)

    n_mcsamples = length(outcome_array)
    @threads for ii in 1:n_mcsamples
        coeff, final_operator, is_valid = montecarlopropagation(circ, pstr, thetas, split_probabilities; is_reversed=true, kwargs...)

        # multiply the coefficient of the backpropagated Pauli with the overlap with the initial state (#TODO: provide initial state)
        # and then multply with (1 - is_valid) to get the final outcome.
        # if not valid, then the error is coeff * overlapwithzero(final_operator), if valid, then the error is 0
        outcome_array[ii] = coeff * overlapwithzero(final_operator) * !is_valid
    end

    return mean(outcome_array)

end

function montecarlopropagation(circ, pstr::PauliString, thetas::AbstractArray, split_probabilities::AbstractArray; is_reversed=false, max_weight=Inf, max_freq=Inf, max_sins=Inf, kwargs...)
    # Reverse the circ if it is not already done. Allocates memory.
    if !is_reversed
        circ = reverse(circ)
    end
    param_idx = length(thetas)
    prob_idx = length(split_probabilities)

    is_valid = true

    oper = pstr.operator
    coeff = pstr.coeff
    for gate in circ
        coeff, oper = mcapply(gate, oper, coeff, thetas[param_idx], split_probabilities[prob_idx]; kwargs...)

        # check truncations
        # TODO: make this properly
        if truncateweight(oper, max_weight)
            is_valid = false
            oper = typeof(oper)(0)
            break
        end
        if truncatefrequency(coeff, max_freq)
            is_valid = false
            oper = typeof(oper)(0)
            break
        end
        if truncatesins(coeff, max_sins)
            is_valid = false
            oper = typeof(oper)(0)
            break
        end

        # TODO: Custom truncation functions?

        # decrement the parameter index if gate is parametrized
        if isa(gate, ParametrizedGate)
            param_idx -= 1
        end
        # always decrement the probability index
        prob_idx -= 1
    end
    return coeff, oper, is_valid
end

# if we don't define an mcapply() function, assume that it doesn't split and that the apply function does the job
mcapply(gate, oper, theta, coeff=1.0, args...) = apply(gate, oper, theta, coeff)

# mcapply for Pauli gates
function mcapply(gate::PauliGateUnion, oper, theta, coeff=1.0, split_prob=0.5)  # 

    if !commutes(gate, oper)
        # if the gate does not commute with the operator, split with probability split_prob int a random direction.
        if rand() > split_prob
            # TODO: make sure this works in general and not just for numerical coefficients. For example, add to freq and nsins.
            sign, oper = getnewoperator(gate, oper)
        end
    end
    # the sign has abs(sign) = 1, so updating the coefficient with abs(sign) doesn't do anything
    return coeff, oper
end

function mcapply(gate::AmplitudeDampingNoise, args...)  # 
    # TODO
    println("AmplitudeDampingNoise not implemented yet")
end
