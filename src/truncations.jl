# Custom truncation function

# Define the custom truncation functions by dissipation-assisted damping
"""
    truncatedampingcoeff(
        pstr::Integer, 
        coeff::Real, 
        gamma::Real, 
        min_abs_coeff::Float64
    )

Custom truncation function with dissipation-assisted damping of coefficients.

This function damps the coefficient `coeff` scaling with the weight of an interger Pauli string `pstr`. 
The damping factor is `gamma`. 
If the coefficient, damped by an exponential factor, falls below `min_abs_coeff`, 
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
function truncatedampingcoeff(
    pstr::PauliStringType, coeff, gamma::Real, min_abs_coeff::Real
)

    return abs(tonumber(coeff)) * 10.0^(-gamma * countweight(pstr)) < min_abs_coeff
end


## TODO: Make actual use of this file or remove.
## from here on, the are just for the backend

## Truncations on the Pauli string:

# Return `true` if an integer Pauli string should be truncated because its weight (i.e., number of non-identity Paulis) is larger than `max_weight`. 
function truncateweight(pstr::PauliStringType, max_weight::Real)
    return countweight(pstr) > max_weight
end


## Truncations on the coefficient:

# Truncations on unsuitable coefficient types defaults to false.
function truncatemincoeff(coeff, min_abs_coeff::Real)
    return false
end


# Return `true` if `abs(coeff) < min_abs_coeff`. 
function truncatemincoeff(coeff::Real, min_abs_coeff::Real)
    return abs(coeff) < min_abs_coeff
end


# Return `true` if `abs(path_property.coeff) < min_abs_coeff`. 
function truncatemincoeff(path_property::PProp, min_abs_coeff::Real) where {PProp<:PathProperties}
    if hasfield(PProp, :coeff)
        return abs(path_property.coeff) < min_abs_coeff
    else
        return false
    end
end


# Return `true` if  `PathProperties.freq > max_freq`. Truncations on coefficients should default to false if it is not applicable for a type.
function truncatefrequency(coeff, max_freq::Real)
    return false
end


# Return `true` if  `path_properties.freq > max_freq`.
function truncatefrequency(path_properties::PProp, max_freq::Real) where {PProp<:PathProperties}
    if hasfield(PProp, :freq)
        return path_properties.freq > max_freq
    else
        return false
    end
end

# Return `true` if  `PathProperties.nsins > max_sins`. Truncations on coefficients should default to false if it is not applicable for a type.
function truncatesins(coeff, max_sins::Real)
    return false
end


# Return `true` if  `path_properties.nsins > max_sins`.
function truncatesins(path_properties::PProp, max_sins::Real) where {PProp<:PathProperties}
    if hasfield(PProp, :nsins)
        return path_properties.nsins > max_sins
    else
        return false
    end
end

