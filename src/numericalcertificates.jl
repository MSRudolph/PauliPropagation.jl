using Statistics



"""
    estimateaverageerror(circ, pstr::PauliString, n_mcsamples::Integer; initialstatefunc=overlapwithzero, circuit_reversed=false, kwargs...)

Function to estimate the average error of a truncated Pauli propagation simulation using Monte Carlo sampling.
This version of the function is only valid for a circuit consisting of `PauliGates` non-parametrized non-splitting gates (for example Cliffords).
Returns the mean squared error of the truncated Pauli propagation simulation averaged over the angles of all `PauliGate`s in the interval [-pi, pi].

An initial state overlap function `initialstatefunc` can be provided to calculate the overlap of the backpropagated Pauli operators with the initial state.
Importantly, the `kwargs` can be used to set the truncation parameters of the Pauli propagation.
"""
function estimateaverageerror(circ, pstr::PauliString, n_mcsamples::Integer; initialstatefunc=overlapwithzero, circuit_reversed=false, kwargs...)
    # this function is only valid for PauliGates and non-splitting non-parametrized gates
    # this will not not error for non-parametrized splitting gates, e.g. T-gates. 
    if !all(g -> isa(g, PauliGateUnion) || isa(g, StaticGate), circ)
        throw("`circ` must contain only (Fast)PauliGates and CliffordGates. Otherwise provide `thetas` manually.")
    end

    # integration range over thetas is [-pi, pi]
    thetas = π

    return estimateaverageerror(circ, pstr, n_mcsamples, thetas; initialstatefunc=initialstatefunc, circuit_reversed=circuit_reversed, kwargs...)
end



"""
    estimateaverageerror(circ, pstr::PauliString, n_mcsamples::Integer, thetas; initialstatefunc=overlapwithzero, circuit_reversed=false, kwargs...)

Function to estimate the average error of a truncated circuit simulation using Monte Carlo sampling.
Returns the mean squared error of the truncated Pauli propagation simulation averaged over the `thetas`∈ [theta, theta] of the angle `theta` of each `PauliGate`.

The length the `thetas` vector should be equal to the number of parametrized gates in the circuit. Alternative, `thetas` can be a single real number applicable for all parametrized gates.
For `PauliGates`, the value should be the integration range of the parameters around zero.
For other currently supported parametrized gates, potential splitting probabilities can be derived from the parameters (e.g., for `AmplitudeDampingNoise`). 
We currently support no non-parametrized splitting gates. This may change in the future.

An initial state overlap function `initialstatefunc` can be provided to calculate the overlap of the backpropagated Pauli operators with the initial state.
Importantly, the `kwargs` can be used to set the truncation parameters of the Pauli propagation.
"""
function estimateaverageerror(circ, pstr::PauliString, n_mcsamples::Integer, thetas; initialstatefunc=overlapwithzero, circuit_reversed=false, kwargs...)
    # length(thetas) should be equal to the number of parametrized gates in the circuit

    split_probabilities = _calculatesplitprobabilities(circ, thetas)

    outcome_array = zeros(n_mcsamples)

    return estimateaverageerror!(circ, pstr, outcome_array, thetas, split_probabilities; initialstatefunc=initialstatefunc, circuit_reversed=circuit_reversed, kwargs...)

end

"""
    estimateaverageerror!(circ, pstr::PauliString, outcome_array::AbstractVector, thetas, split_probabilities; initialstatefunc=paulioverlapwithzero, circuit_reversed=false, kwargs...)

In-place version of `estimateaverageerror`. This function takes an array of length n_mcsamples as an argument and modifies it in-place. 
It further assumes that the `thetas` and `split_probabilities` are already correctly calculated and provided as arguments. 
In general they will be vectors, but they can also be real numbers.
"""
function estimateaverageerror!(circ, pstr::PauliString, outcome_array::AbstractVector, thetas, split_probabilities; initialstatefunc=paulioverlapwithzero, circuit_reversed=false, kwargs...)
    # This function takes an outcome_array as an argument and modifies it in-place.

    # length(thetas) should be equal to the number of parametrized gates in the circuit
    if isa(thetas, AbstractVector)
        if length(thetas) != countparameters(circ)
            throw("Vector `thetas` must have same length the number of parametrized gates in `circ`.")
        end
    end
    # length(split_probabilities) should be equal to the number of total gates in the circuit
    if isa(split_probabilities, AbstractArray)
        if length(split_probabilities) != length(circ)
            throw("Vector `split_probabilities` must have same length as circuit.")
        end
    end

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
        outcome_array[ii] = getnumcoeff(final_pstr.coeff) * initialstatefunc(final_pstr.operator) * (1 - is_valid)
    end

    return mean(outcome_array)

end



"""
    montecarlopropagation(circ, pstr::PauliString; circuit_reversed=false, max_weight=Inf, max_freq=Inf, max_sins=Inf, kwargs...)

Perform a single Monte Carlo propagation of a Pauli string through a circuit.
Assumes that the circuit contains only (Fast)PauliGates and CliffordGates; and sets the integration rage to [-pi, pi].
"""
function montecarlopropagation(circ, pstr::PauliString; circuit_reversed=false, max_weight=Inf, max_freq=Inf, max_sins=Inf, kwargs...)
    if !all(g -> isa(g, PauliGateUnion) || isa(g, StaticGate), circ)
        throw("`circ` must contain only (Fast)PauliGates and CliffordGates. Otherwise provide `thetas` manually.")
    end

    thetas = π

    return montecarlopropagation(
        circ, pstr, thetas;
        circuit_reversed=circuit_reversed, max_weight=max_weight, max_freq=max_freq, max_sins=max_sins, kwargs...
    )
end

"""
    montecarlopropagation(circ, pstr::PauliString, thetas; circuit_reversed=false, max_weight=Inf, max_freq=Inf, max_sins=Inf, kwargs...)

Perform a single Monte Carlo propagation of a Pauli string through a circuit.

The length the `thetas` vector should be equal to the number of parametrized gates in the circuit. Alternative, `thetas` can be a single real number applicable for all parametrized gates.
For `PauliGates`, the value should be the integration range of the parameters around zero.
For other currently supported parametrized gates, potential splitting probabilities can be derived from the parameters (e.g., for `AmplitudeDampingNoise`). 
We currently support no non-parametrized splitting gates. This may change in the future.
"""
function montecarlopropagation(circ, pstr::PauliString, thetas; circuit_reversed=false, max_weight=Inf, max_freq=Inf, max_sins=Inf, kwargs...)
    split_probabilities = _calculatesplitprobabilities(circ, thetas)
    return montecarlopropagation(
        circ, pstr, thetas, split_probabilities;
        circuit_reversed=circuit_reversed, max_weight=max_weight, max_freq=max_freq, max_sins=max_sins, kwargs...
    )
end

"""
    montecarlopropagation(circ, pstr::PauliString, thetas, split_probabilities; circuit_reversed=false, max_weight=Inf, max_freq=Inf, max_sins=Inf, kwargs...)

Perform a single Monte Carlo propagation of a Pauli string through a circuit.

It further assumes that the `thetas` and `split_probabilities` are already correctly calculated and provided as arguments. 
In general they will be vectors, but they can also be real numbers.
"""
function montecarlopropagation(circ, pstr::PauliString, thetas, split_probabilities; circuit_reversed=false, max_weight=Inf, max_freq=Inf, max_sins=Inf, kwargs...)
    # Reverse the circ if it is not already done. Allocates memory.
    if circuit_reversed
        circ = reverse(circ)
    end
    param_idx = length(thetas)
    prob_idx = length(split_probabilities)

    is_valid = true

    pauli = pstr.operator
    coeff = copy(pstr.coeff)
    for gate in circ

        coeff, pauli = mcapply(gate, pauli, coeff, _getelmt(thetas, param_idx), _getelmt(split_probabilities, prob_idx); kwargs...) # TODO: this currently allocates a lot of memory.

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

## Monte Carlo apply functions

"""
    mcapply(gate, pauli, coeff, theta, split_probability; kwargs...)

Default `mcapply` function for numerical certificates. Assumes that the `apply` function for `gate` does the job.
"""
mcapply(gate, pauli, coeff, theta, split_probability; kwargs...) = apply(gate, pauli, theta, coeff; kwargs...) # TODO: have more sensible default functions once the numerical certificate is figured out.

"""
    mcapply(gate::PauliGateUnion, pauli, coeff, theta, split_prob=0.5; kwargs...) 

MC apply function for Pauli gates. If the gate commutes with the operator, the operator is left unchanged. Else the operator is split off with a probability 1 - `split_prob`.
"""
function mcapply(gate::PauliGateUnion, pauli, coeff, theta, split_prob=0.5; kwargs...)

    if !commutes(gate, pauli)
        # if the gate does not commute with the operator, remain with probability `split_prob` and split off with probability 1 - `split_prob`.
        if rand() > split_prob
            # branch into the new Pauli
            sign, pauli = getnewoperator(gate, pauli)
            # for PathProperties: increment sin and frequency count
            _incrementsinandfreq!(coeff)
        else
            # Pauli doesn't get changed
            # for PathProperties: increment cos and frequency count
            _incrementcosandfreq!(coeff)
        end
    end
    # the sign has abs(sign) = 1, so updating the coefficient with abs(sign) doesn't do anything
    return coeff, pauli
end

"""
TODO: `mcapply` function for AmplitudeDampingNoise.
"""
function mcapply(gate::AmplitudeDampingNoise, args...; kwargs...)  # 
    # TODO
    throw("AmplitudeDampingNoise not implemented yet")
end



### Utilities
_getelmt(arr::AbstractArray, idx::Integer) = arr[idx]
_getelmt(num::Number, idx::Integer) = num

"""
Function that automatically calculates the vector of splitting probabilities of the gates in the circuit based on a vector of thetas.
For Pauli gates, the theta value is interpreted as the limits of the integration [-theta, theta].
For AmplitudeDampingNoise, the splitting probability is the damping rate.
"""
function _calculatesplitprobabilities(circ::AbstractArray, thetas::AbstractArray)
    if length(thetas) != countparameters(circ)
        throw("Vector `thetas` must have same length the number of parametrized gates in `circ`.")
    end

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

"""
Function that automatically calculates the splitting probability of the gates in the circuit based on a one number theta.
For Pauli gates, the theta value is interpreted as the limits of the integration [-theta, theta].
For AmplitudeDampingNoise, the splitting probability is the damping rate.
"""
function _calculatesplitprobabilities(circ::AbstractArray, theta::Number)
    return 0.5 * (1 + sin(2 * theta) / (2 * theta))
end