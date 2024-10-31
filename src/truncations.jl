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
function truncatesins(path_properties, max_sins::Real)
    return false
end

function truncatesins(path_properties::PathProperties, max_sins::Real)
    return path_properties.nsins > max_sins
end

# Custom truncation function

# Define the custom truncation functions by dissipation-assisted damping
function truncatedampingcoeff(
    pstr::PauliStringType, coeff::Float64, gamma::Float64, min_abs_coeff
)
  return abs(coeff) * exp(- gamma * countweight(pstr))  < min_abs_coeff
end

function customtruncatedampingcoeff(
    pstr_dict, gamma::Float64, min_abs_coeff
)
"""
Custom truncation function with dissipation-assisted damping of coefficients.

Truncate Pauli strings with coefficients < `min_abs_coeff` and damping `gamma`:
c(P) * exp(- gamma * w(P)) < `min_abs_coeff`

Args:
    pstr_dict (dict): Dictionary of Pauli strings.
    gamma (float): Damping gamma.
    min_abs_coeff (float): Minimal absolute value of the coefficient.

Returns:
    None: Changes `pstr_dict` in-place.
"""
  for (pstr, coeff) in pstr_dict
      if truncatedampingcoeff(pstr, coeff, gamma, min_abs_coeff)
          delete!(pstr_dict, pstr)
      end
  end
  return
end
