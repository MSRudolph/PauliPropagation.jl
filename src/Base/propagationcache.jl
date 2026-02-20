###
##
# Propagation caches carry the main term sum and auxiliary data structures for efficient propagation.
##
###

abstract type AbstractPropagationCache end

## The interface to hook into with each library
PropagationCache(thing::AbstractTermSum) = _thrownotimplemented(thing, :PropagationCache)

mainsum(prop_cache::AbstractPropagationCache) = _thrownotimplemented(prop_cache, :mainsum)
auxsum(prop_cache::AbstractPropagationCache) = _thrownotimplemented(prop_cache, :auxsum)

setmainsum!(prop_cache::AbstractPropagationCache, new_mainsum) = _thrownotimplemented(prop_cache, :setmainsum!)
setauxsum!(prop_cache::AbstractPropagationCache, new_auxsum) = _thrownotimplemented(prop_cache, :setauxsum!)

function activesum(prop_cache::AbstractPropagationCache)
    _thrownotimplemented(prop_cache, :activesum)
end

function swapsums!(prop_cache::AbstractPropagationCache)
    temp = mainsum(prop_cache)
    setmainsum!(prop_cache, auxsum(prop_cache))
    setauxsum!(prop_cache, temp)
    return prop_cache
end

function copyswapsums!(prop_cache::AbstractPropagationCache)
    # instead of just swapping auxsum(prop_cache) and mainsum(prop_cache),
    # we need to copy the mainsum into the aux sum and then swap the sums
    copy!(auxsum(prop_cache), mainsum(prop_cache))
    swapsums!(prop_cache)
    return prop_cache
end

StorageType(prop_cache::AbstractPropagationCache) = StorageType(mainsum(prop_cache))

nsites(prop_cache::AbstractPropagationCache) = nsites(mainsum(prop_cache))

function Base.length(prop_cache::AbstractPropagationCache)
    return _length(StorageType(prop_cache), prop_cache)
end
function _length(::DictStorage, prop_cache::AbstractPropagationCache)
    return length(mainsum(prop_cache))
end

function _length(::ArrayStorage, prop_cache::AbstractPropagationCache)
    return activesize(prop_cache)
end

Base.isempty(prop_cache::AbstractPropagationCache) = length(prop_cache) == 0

function capacity(prop_cache::AbstractPropagationCache)
    return _capacity(StorageType(prop_cache), prop_cache)
end

function _capacity(::DictStorage, prop_cache::AbstractPropagationCache)
    return length(storage(mainsum(prop_cache)).slots)
end

function _capacity(::ArrayStorage, prop_cache::AbstractPropagationCache)
    return length(mainsum(prop_cache))
end

## Interface for vector-based propagation caches

activesize(prop_cache::AbstractPropagationCache) = _thrownotimplemented(prop_cache, :activesize)
setactivesize!(prop_cache::AbstractPropagationCache, new_size::Int) = _thrownotimplemented(prop_cache, :setactivesize!)


activeterms(prop_cache::AbstractPropagationCache) = view(terms(mainsum(prop_cache)), 1:activesize(prop_cache))
activecoeffs(prop_cache::AbstractPropagationCache) = view(coefficients(mainsum(prop_cache)), 1:activesize(prop_cache))
activeauxterms(prop_cache::AbstractPropagationCache) = view(terms(auxsum(prop_cache)), 1:activesize(prop_cache))
activeauxcoeffs(prop_cache::AbstractPropagationCache) = view(coefficients(auxsum(prop_cache)), 1:activesize(prop_cache))

flags(prop_cache::AbstractPropagationCache) = _thrownotimplemented(prop_cache, :flags)
indices(prop_cache::AbstractPropagationCache) = _thrownotimplemented(prop_cache, :indices)
activeflags(prop_cache::AbstractPropagationCache) = view(flags(prop_cache), 1:activesize(prop_cache))
activeindices(prop_cache::AbstractPropagationCache) = view(indices(prop_cache), 1:activesize(prop_cache))
lastactiveindex(prop_cache::AbstractPropagationCache) = activeindices(prop_cache)[end]



function Base.resize!(prop_cache::AbstractPropagationCache, new_size::Int)
    _thrownotimplemented(prop_cache, :resize!)
end

## Back-conversions 
function (::Type{TS})(prop_cache::AbstractPropagationCache) where TS<:AbstractTermSum
    return convert(TS, _cachetosum!(StorageType(prop_cache), prop_cache))
end

_cachetosum!(::DictStorage, prop_cache::AbstractPropagationCache) = mainsum(prop_cache)

function _cachetosum!(::ArrayStorage, prop_cache::AbstractPropagationCache)
    # convert back to TermSum
    term_sum = mainsum(prop_cache)
    term_sum = resize!(term_sum, activesize(prop_cache))
    return term_sum
end

function _cachetosum!(::Type{ST}, ::PC) where {ST<:StorageType,PC<:AbstractPropagationCache}
    throw(ErrorException("cachetosum!(::$(ST), ::$(PC)) not implemented."))
end
