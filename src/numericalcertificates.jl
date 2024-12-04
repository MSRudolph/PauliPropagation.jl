using Statistics

"""
    estimateaverageerror(circ, pstr::PauliString, n_mcsamples::Integer, thetas=π; stateoverlapfunc=overlapwithzero, circuit_reversed=false, kwargs...)

Function to estimate the average error of a truncated circuit simulation using Monte Carlo sampling.
Returns the mean squared error of the truncated Pauli propagation simulation averaged over the `thetas`∈ [theta, theta] of the angle `theta` of each `PauliGate`.

The length the `thetas` vector should be equal to the number of parametrized gates in the circuit. 
Alternatively, `thetas` can be a single real number applicable for all parametrized gates.
The default `thetas=π` or any other non-array values assume that the circuit consists only of `(Fast)PauliGates` -`CliffordGates`.
For `PauliGates`, the value should be the integration range of the parameters around zero.
For other currently supported parametrized gates, potential splitting probabilities can be derived from the parameters (e.g., for `AmplitudeDampingNoise`). 
We currently support no non-parametrized splitting gates. This may change in the future.

An initial state overlap function `stateoverlapfunc` can be provided to calculate the overlap of the backpropagated Pauli strings with the initial state.
Importantly, the `kwargs` can be used to set the truncation parameters of the Pauli propagation. Currently supported are `max_weight`, `max_freq`, and `max_sins`.
Note that `min_abs_coeff` is not supported here, as we consider errors integrated over the angles. `max_freq` effectively truncates small coefficients below (1/2)^`max_freq` on average over `thetas ∈ [-π, π]`.
"""
function estimateaverageerror(circ, pstr::PauliString, n_mcsamples::Integer, thetas=π; stateoverlapfunc=overlapwithzero, circuit_reversed=false, kwargs...)
    # this function is only valid for ParametrizedGates and non-splitting non-parametrized gates (here only CliffordGates).
    # this will not not error for non-parametrized splitting gates, e.g. T-gates. 
    if !all(g -> isa(g, PauliGateUnion) || isa(g, CliffordGate), circ)
        throw("`circ` must contain only `ParameterizedGate`s and `CliffordGate`s.")
    end

    split_probabilities = _calculatesplitprobabilities(circ, thetas)

    error_array = zeros(n_mcsamples)

    return estimateaverageerror!(circ, pstr, error_array, thetas, split_probabilities; stateoverlapfunc=stateoverlapfunc, circuit_reversed=circuit_reversed, kwargs...)

end

"""
    estimateaverageerror!(circ, pstr::PauliString, error_array::AbstractVector, thetas, split_probabilities; stateoverlapfunc=overlapwithzero, circuit_reversed=false, kwargs...)

In-place version of `estimateaverageerror`. This function takes an array `error_array` of length `n_mcsamples` as an argument and modifies it in-place. 
It further assumes that the `thetas` and `split_probabilities` are already correctly calculated and provided as arguments. 
In general they will be vectors, but they can also be real numbers.
"""
function estimateaverageerror!(circ, pstr::PauliString, error_array::AbstractVector, thetas, split_probabilities; stateoverlapfunc=overlapwithzero, circuit_reversed=false, kwargs...)
    # This function takes an error_array as an argument and modifies it in-place.

    # length(thetas) should be equal to the number of parametrized gates in the circuit
    if isa(thetas, AbstractVector)
        if length(thetas) != countparameters(circ)
            throw("Vector `thetas` of length $(length(thetas)) must have same length the number of parametrized gates $(countparameters(circ)) in `circ`.")
        end
    end
    # length(split_probabilities) should be equal to the number of total gates in the circuit
    if isa(split_probabilities, AbstractArray)
        if length(split_probabilities) != length(circ)
            throw("Vector `split_probabilities` of length $(length(split_probabilities)) must have same length as `circ` of length $(length(circ)).")
        end
    end

    # reverse the circuit once 
    if circuit_reversed
        circ = reverse(circ)
    end

    n_mcsamples = length(error_array)
    @threads for ii in 1:n_mcsamples
        final_pstr, is_truncated = montecarlopropagation(circ, pstr, thetas, split_probabilities; circuit_reversed=true, kwargs...)

        # multiply the coefficient of the backpropagated Pauli with the overlap with the initial state
        # and then multply with `is_truncated` to get the final error.
        # if truncated, then the error is coeff * initial_state_func(final_operator), if not truncated, then the error is 0
        error_array[ii] = getnumcoeff(final_pstr.coeff)^2 * stateoverlapfunc(final_pstr.operator)^2 * is_truncated
    end

    return mean(error_array)

end


"""
    montecarlopropagation(circ, pstr::PauliString, thetas=π; circuit_reversed=false, max_weight=Inf, max_freq=Inf, max_sins=Inf, kwargs...)

Perform a single Monte Carlo propagation of a Pauli string through a circuit. Returns the final Pauli string and a boolean indicating whether the path was truncated.

The length the `thetas` vector should be equal to the number of parametrized gates in the circuit. Alternatively, `thetas` can be a single real number applicable for all parametrized gates.
The default `thetas=π` or any other non-array values assume that the circuit consists only of `(Fast)PauliGates` -`CliffordGates`.
For `PauliGates`, `thetas` should be the integration range of the parameters around zero.
For other currently supported parametrized gates, potential splitting probabilities can be derived from the parameters (e.g., for `AmplitudeDampingNoise`). 
We currently support no non-parametrized splitting gates. This may change in the future.
"""
function montecarlopropagation(circ, pstr::PauliString, thetas=π; circuit_reversed=false, max_weight=Inf, max_freq=Inf, max_sins=Inf, kwargs...)
    if !all(g -> isa(g, PauliGateUnion) || isa(g, CliffordGate), circ)
        throw("`circ` must contain only `ParametrizedGate`s and CliffordGates.")
    end

    split_probabilities = _calculatesplitprobabilities(circ, thetas)
    return montecarlopropagation(
        circ, pstr, thetas, split_probabilities;
        circuit_reversed=circuit_reversed, max_weight=max_weight, max_freq=max_freq, max_sins=max_sins, kwargs...
    )
end

"""
    montecarlopropagation(circ, pstr::PauliString, thetas, split_probabilities; circuit_reversed=false, max_weight=Inf, max_freq=Inf, max_sins=Inf, kwargs...)

Perform a single Monte Carlo propagation of a Pauli string through a circuit. Returns the final Pauli string and a boolean indicating whether the path was truncated.

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

    is_truncated = false

    pauli = pstr.operator
    coeff = copy(pstr.coeff)
    for gate in circ

        pauli, coeff = mcapply(gate, pauli, coeff, _getelmt(thetas, param_idx), _getelmt(split_probabilities, prob_idx); kwargs...) # TODO: this currently allocates a lot of memory.

        # check truncations
        # TODO: make this properly
        if truncateweight(pauli, max_weight)
            is_truncated = true
        elseif truncatefrequency(coeff, max_freq)
            is_truncated = true
        elseif truncatesins(coeff, max_sins)
            is_truncated = true
        end

        # TODO: Custom truncation functions?

        # decrement the parameter index if gate is parametrized
        if isa(gate, ParametrizedGate)
            param_idx -= 1
        end
        # always decrement the probability index
        prob_idx -= 1
    end
    return PauliString{typeof(pstr).parameters...}(pstr.nqubits, pauli, coeff), is_truncated
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
            pauli, sign = getnewpaulistring(gate, pauli) # TODO: This allocates
            # for PathProperties: increment sin and frequency count
            coeff = _incrementsinandfreq(coeff)
        else
            # Pauli doesn't get changed
            # for PathProperties: increment cos and frequency count
            coeff = _incrementcosandfreq(coeff)
        end
    end
    # the sign has abs(sign) = 1, so updating the coefficient with abs(sign) doesn't do anything
    return pauli, coeff
end

"""
TODO: `mcapply` function for AmplitudeDampingNoise.
"""
function mcapply(gate::AmplitudeDampingNoise, args...; kwargs...)  # 
    # TODO
    throw("AmplitudeDampingNoise not implemented yet")
end

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
        if isa(gate, ParametrizedGate)
            split_probabilities[ii] = _calculatesplitprobabilities(gate, thetas[theta_idx])
            theta_idx += 1
        end
    end
    return split_probabilities
end

"""
Function that automatically calculates the splitting probability of the gates in the circuit based on a one number theta.
This assumes that the circuit consists only of `(Fast)PauliGates` -`CliffordGates`.
"""
function _calculatesplitprobabilities(circ::AbstractArray, r::Number)
    return 0.5 * (1 + sin(2 * r) / (2 * r))
end

_calculatesplitprobabilities(gate::PauliGateUnion, theta::Number) = 0.5 * (1 + sin(2theta) / (2theta))

_calculatesplitprobabilities(gate::AmplitudeDampingNoise, theta::Number) = theta

_calculatesplitprobabilities(gate, theta::Number) = 0.0


### Utilities
_getelmt(arr::AbstractArray, idx::Integer) = idx > 0 ? arr[idx] : eltype(arr)(0)
_getelmt(num::Number, idx::Integer) = num
