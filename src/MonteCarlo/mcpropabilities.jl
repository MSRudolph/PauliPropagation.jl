### mcprobabilities.jl
##
# This file contains the functions for computing MC probabilities. 
##
###

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
