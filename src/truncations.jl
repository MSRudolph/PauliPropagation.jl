"""
    Truncation function should return 'true' if a path should continue to propagate
    and 'false' if it should be truncated.
"""

function truncateweight(oper, max_weight)
    return countweight(oper) > max_weight
end



"""
Return 'true' if abs(coeff) < min_abs_coeff
"""
function truncatemincoeff(coeff, min_abs_coeff)
    return false
end

function truncatemincoeff(coeff::Float64, min_abs_coeff)
    return abs(coeff) < min_abs_coeff
end

function truncatemincoeff(node::NumericPathProperties, min_abs_coeff)
    return abs(node.coeff) < min_abs_coeff
end


"""
Return 'true' if freq > max_freq
"""
function truncatefrequency(coeff, max_freq::Real)
    return false
end

function truncatefrequency(path_properties::T, max_freq::Real) where {T<:PathProperties}
    return path_properties.freq > max_freq
end



"""
Return 'true' if n_sins > max_sins
"""
function truncatesins(coeff, max_sins::Real)
    return false
end

function truncatesins(path_properties::PathProperties, max_sins::Real)
    return path_properties.nsins > max_sins
end

# Custom truncation function

# Define the custom truncation functions by dissipation-assisted damping

function truncatedampingcoeff(
    pstr::PauliStringType, coeff::Float64, gamma::Float64, min_abs_coeff
)::Bool
"""
    truncatedampingcoeff(
        pstr::PauliStringType, 
        coeff::Float64, 
        gamma::Float64, 
        min_abs_coeff::Float64
    ) -> Bool

Custom truncation function with dissipation-assisted damping of coefficients.

This function damps the coefficient by the weight of a Pauli string `pstr`. 
Its associated `coeff`, and a damping factor `gamma`.If the coefficient, 
damped by an exponential factor, falls below `min_abs_coeff`, 
the function returns `true` to indicate truncation.

# Arguments
- `pstr::PauliStringType`: The Pauli string whose coefficient may be damped.
- `coeff::Float64`: The coefficient associated with the Pauli string `pstr`.
- `gamma::Float64`: Gamma, rate of exponential decay in the damping process.
- `min_abs_coeff::Float64`: The minimum value of the coefficient for truncation.

# Returns
- `Bool`: `true` if the damped coefficient of `pstr` < `min_abs_coeff`;
  `false` otherwise.

# Details
The function evaluates the condition:
`abs(coeff) * exp(-gamma * w(pstr)) < min_abs_coeff`

where `w(pstr)` is the weight of the Pauli string (computed by `countweight`). 
The damping factor `gamma` controls the exponential decay.

# Examples
```julia
truncatedampingcoeff(pstr, 0.8, 0.5, 0.01)
"""
  return abs(coeff) * exp(- gamma * countweight(pstr))  < min_abs_coeff
end
