## TODO: Make actual use of this ile or remove.

"""
Return `true` if a Pauli string in its integer representation should be truncated because its weight (i.e., number of non-identity Paulis) is larger than `max_weight`. 
"""
function truncateweight(oper::PauliStringType, max_weight::Real)
    return countweight(oper) > max_weight
end

"""
Return `true` if `abs(coeff) < min_abs_coeff`. Truncations on coefficients should default to false if it is not applicable for a type.
"""
function truncatemincoeff(coeff, min_abs_coeff::Real)
    return false
end

function truncatemincoeff(coeff::Float64, min_abs_coeff::Real)
    return abs(coeff) < min_abs_coeff
end

function truncatemincoeff(node::NumericPathProperties, min_abs_coeff::Real)
    return abs(node.coeff) < min_abs_coeff
end


"""
Return `true` if  `PathProperties.freq > max_freq`. Truncations on coefficients should default to false if it is not applicable for a type.
"""
function truncatefrequency(coeff, max_freq::Real)
    return false
end

function truncatefrequency(path_properties::T, max_freq::Real) where {T<:PathProperties}
    return path_properties.freq > max_freq
end



"""
Return `true` if  `PathProperties.nsins > max_sins`. Truncations on coefficients should default to false if it is not applicable for a type.
"""
function truncatesins(path_properties, max_sins::Real)
    return false
end

function truncatesins(path_properties::PathProperties, max_sins::Real)
    return path_properties.nsins > max_sins
end

