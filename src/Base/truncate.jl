function truncate!(term_sum::AbstractTermSum; min_abs_coeff::Real, customtruncationfunc::F=_alwaysfalse) where F<:Function
    prop_cache = truncate!(truncfunc, PropagationCache(term_sum))
    return mainsum(prop_cache)
end

function truncate!(prop_cache::AbstractPropagationCache; min_abs_coeff::Real=0.0, customtruncationfunc::F=_alwaysfalse) where F<:Function

    # bundle truncation functions 
    truncfunc = (pstr, coeff) -> truncatemincoeff(coeff, min_abs_coeff) || customtruncationfunc(pstr, coeff)

    prop_cache = truncate!(truncfunc, prop_cache)
    return prop_cache
end

# this can can be a term sum or a propagation cache
function truncate!(truncfunc::F, prop_object) where F<:Function
    return truncate!(StorageType(prop_object), truncfunc, prop_object)
end


function truncate!(::DictStorage, truncfunc::F, prop_cache::AbstractPropagationCache; kwargs...) where F<:Function
    term_sum = mainsum(prop_cache)
    term_sum = truncate!(StorageType(term_sum), truncfunc, term_sum; kwargs...)
    setmainsum!(prop_cache, term_sum)
    return prop_cache
end

function truncate!(::DictStorage, truncfunc::F, term_sum::AbstractTermSum) where F<:Function
    for (pstr, coeff) in term_sum
        if truncfunc(pstr, coeff)
            delete!(term_sum, pstr)
        end
    end
    return term_sum
end

function truncate!(::ArrayStorage, truncfunc::F, prop_cache::AbstractPropagationCache) where F<:Function

    if isempty(prop_cache)
        return prop_cache
    end

    # TODO: view and active size interface for vector sums

    # flag the indices that we keep
    keepfunc(pstr, coeff) = !truncfunc(pstr, coeff)
    flag!(keepfunc, prop_cache)

    filterviaflags!(prop_cache)

    return prop_cache
end

function truncate!(::ArrayStorage, truncfunc::F, term_sum::AbstractTermSum) where F<:Function
    # convert to propagation cache for easier handling
    prop_cache = PropagationCache(term_sum)

    prop_cache = truncate!(ArrayStorage(), truncfunc, prop_cache)

    return mainsum(prop_cache)
end


# Truncations on unsuitable coefficient types defaults to false.
function truncatemincoeff(coeff, min_abs_coeff)
    return false
end


# This should work for any complex and real coefficient
function truncatemincoeff(coeff::Number, min_abs_coeff::Real)
    return abs(coeff) < min_abs_coeff
end


_alwaysfalse(::Any...) = false
