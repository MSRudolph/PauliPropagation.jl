
function PauliPropagation.overlapwithzero(prop_cache::PropagationCache)
    total = zero(coefftype(prop_cache))
    for ii in eachindex(viewterms(prop_cache))
        total += tonumber(prop_cache.coeffs[ii]) * !containsXorY(prop_cache.terms[ii])
    end
    return total
end


function PauliPropagation.overlapwithcomputational(prop_cache::PropagationCache, onebitinds)
    val = zero(coefftype(prop_cache))
    for i in eachindex(viewterms(prop_cache))
        val += tonumber(prop_cache.coeffs[i]) * PauliPropagation._calcsignwithones(prop_cache.terms[i], onebitinds)
    end
    return val
end