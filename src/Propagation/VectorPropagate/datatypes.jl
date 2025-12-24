###
##
# A PropagationCache carries two VectorPauliSums
# and flags and indices for propagation.
##
###

mutable struct PropagationCache{VT,VC,VB,VI}
    vecpsum::VectorPauliSum{VT,VC}
    aux_vecpsum::VectorPauliSum{VT,VC}
    flags::VB
    indices::VI

    # we will over-allocate the arrays and keep track of the non-empty size
    active_size::Int

    function PropagationCache(
        vecpsum::VectorPauliSum{VT,VC}, aux_vecpsum::VectorPauliSum{VT,VC}, flags::VB, indices::VI, active_size
    ) where {VT,VC,VB,VI}
        #assert that all terms are at least as long as active_size
        @assert length(vecpsum.terms) >= active_size "Length of terms array is less than active_size."
        @assert length(vecpsum.coeffs) >= active_size "Length of coeffs array is less than active_size."
        @assert length(aux_vecpsum.terms) >= active_size "Length of aux_terms array is less than active_size."
        @assert length(aux_vecpsum.coeffs) >= active_size "Length of aux_coeffs array is less than active_size."
        @assert length(flags) >= active_size "Length of flags array is less than active_size."
        @assert length(indices) >= active_size "Length of indices array is less than active_size."

        return new{VT,VC,VB,VI}(vecpsum, aux_vecpsum, flags, indices, active_size)
    end
end

function PropagationCache(vecpsum::VectorPauliSum{VT,VC}) where {VT,VC}
    aux_vecpsum = similar(vecpsum)
    flags = similar(vecpsum.terms, Bool)
    indices = similar(vecpsum.terms, Int)
    return PropagationCache(vecpsum, aux_vecpsum, flags, indices, length(vecpsum.terms))
end

PropagationCache(pstr::PauliString) = PropagationCache(VectorPauliSum(pstr))

function PropagationCache(psum::PauliSum)
    return PropagationCache(VectorPauliSum(psum))
end

function Base.show(io::IO, prop_cache::PropagationCache)
    println(io, "PropagationCache with $(prop_cache.active_size) terms:")
    for i in 1:prop_cache.active_size
        if i > 20
            println(io, "  ...")
            break
        end
        println(io, prop_cache.coeffs[i], " * $(bits(prop_cache.terms[i]))")
    end
end


function Base.resize!(prop_cache::PropagationCache, n_new::Int)
    resize!(prop_cache.vecpsum, n_new)
    resize!(prop_cache.aux_vecpsum, n_new)
    resize!(prop_cache.flags, n_new)
    resize!(prop_cache.indices, n_new)
    return prop_cache
end

termtype(prop_cache::PropagationCache{VT,VC,VB,VI}) where {VT,VC,VB,VI} = eltype(VT)
PauliPropagation.paulitype(prop_cache::PropagationCache) = termtype(prop_cache)
PauliPropagation.coefftype(prop_cache::PropagationCache{VT,VC,VB,VI}) where {VT,VC,VB,VI} = eltype(VC)

viewterms(prop_cache::PropagationCache) = view(prop_cache.vecpsum.terms, 1:prop_cache.active_size)
viewcoeffs(prop_cache::PropagationCache) = view(prop_cache.vecpsum.coeffs, 1:prop_cache.active_size)
viewauxterms(prop_cache::PropagationCache) = view(prop_cache.aux_vecpsum.terms, 1:prop_cache.active_size)
viewauxcoeffs(prop_cache::PropagationCache) = view(prop_cache.aux_vecpsum.coeffs, 1:prop_cache.active_size)
viewflags(prop_cache::PropagationCache) = view(prop_cache.flags, 1:prop_cache.active_size)
viewindices(prop_cache::PropagationCache) = view(prop_cache.indices, 1:prop_cache.active_size)

term(trm::Integer) = trm
term(pstr::PauliString) = term(pstr.term)


Base.length(prop_cache::PropagationCache) = length(prop_cache.vecpsum)
Base.isempty(prop_cache::PropagationCache) = prop_cache.active_size == 0