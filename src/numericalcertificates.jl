function montecarlosampling(circ, oper::Integer, thetas::AbstractArray, split_probabilities::AbstractArray; is_reversed=false, kwargs...)
    # length(thetas) should be equal to the number of parametrized gates in the circuit
    # length(split_probabilities) should be equal to the number of total gates in the circuit

    # Reverse the circ if it is not already done. Allocates memory.
    if !is_reversed
        circ = reverse(circ)
    end
    param_idx = length(thetas)
    prob_idx = length(split_probabilities)

    is_valid = true
    coeff = 1.0
    for gate in circ
        coeff, oper = mcapply(gate, oper, coeff, thetas[param_idx], split_probabilities[prob_idx]; kwargs...)

        # TODO: Check if the operator is truncated.

        # decrement the parameter index if gate is parametrized
        if isa(gate, ParametrizedGate)
            param_idx -= 1
        end
        # always decrement the probability index
        prob_idx -= 1
    end
    return oper, is_valid
end

# if we don't define an mcapply() function, assume that it doesn't split and that the apply function does the job
mcapply(gate, oper, theta, coeff=1.0, args...) = apply(gate, oper, theta, coeff)

function mcapply(gate::PauliGateUnion, oper, theta, coeff=1.0, split_prob=0.5)  # 

    if !commutes(gate, oper)
        # if the gate does not commute with the operator, split with probability split_prob int a random direction.
        if rand() > split_prob
            sign, oper = getnewoperator(gate, oper)
        end
    end
    # the sign has abs(sign) = 1, so updating the coefficient with abs(sign) doesn't do anything
    return coeff, oper
end

# function mcapply(gate::AmplitudeDampingNoise, oper, theta, coeff=1.0, split_prob=0.5)  # 
#     # TODO
# end
