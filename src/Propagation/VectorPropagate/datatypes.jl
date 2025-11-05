mutable struct PropagationCache{VT,VC,VB,VI}
    terms::VT
    coeffs::VC
    aux_terms::VT
    aux_coeffs::VC
    flags::VB
    indices::VI

    # we will over-allocate the arrays and keep track of the non-empty size
    active_size::Int

    function PropagationCache(terms, coeffs, aux_terms, aux_coeffs, flags, indices, active_size)
        #assert that all terms are at least as long as active_size
        @assert length(terms) >= active_size "Length of terms array is less than active_size."
        @assert length(coeffs) >= active_size "Length of coeffs array is less than active_size."
        @assert length(aux_terms) >= active_size "Length of aux_terms array is less than active_size."
        @assert length(aux_coeffs) >= active_size "Length of aux_coeffs array is less than active_size."
        @assert length(flags) >= active_size "Length of flags array is less than active_size."
        @assert length(indices) >= active_size "Length of indices array is less than active_size."

        return new{typeof(terms),typeof(coeffs),typeof(flags),typeof(indices)}(terms, coeffs, aux_terms, aux_coeffs, flags, indices, active_size)
    end
end

function PropagationCache(terms, coeffs)
    aux_terms = similar(terms)
    aux_coeffs = similar(coeffs)
    flags = similar(terms, Bool)
    indices = similar(terms, Int)
    return PropagationCache(terms, coeffs, aux_terms, aux_coeffs, flags, indices, length(terms))
end

PropagationCache(pstr::PauliString) = PropagationCache(PauliSum(pstr))

function PropagationCache(psum::PauliSum)
    nterms = length(psum.terms)
    TT = paulitype(psum)
    CT = coefftype(psum)

    terms = Vector{TT}(undef, nterms)
    coeffs = Vector{CT}(undef, nterms)

    for (i, (trm, cff)) in enumerate(psum.terms)
        terms[i] = trm
        coeffs[i] = cff
    end

    return PropagationCache(terms, coeffs)
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
    resize!(prop_cache.terms, n_new)
    resize!(prop_cache.coeffs, n_new)
    resize!(prop_cache.aux_terms, n_new)
    resize!(prop_cache.aux_coeffs, n_new)
    resize!(prop_cache.flags, n_new)
    resize!(prop_cache.indices, n_new)
    return prop_cache
end

termtype(prop_cache::PropagationCache{VT,VC,VB,VI}) where {VT,VC,VB,VI} = eltype(VT)
PauliPropagation.paulitype(prop_cache::PropagationCache) = termtype(prop_cache)
PauliPropagation.coefftype(prop_cache::PropagationCache{VT,VC,VB,VI}) where {VT,VC,VB,VI} = eltype(VC)

viewterms(prop_cache::PropagationCache) = view(prop_cache.terms, 1:prop_cache.active_size)
viewcoeffs(prop_cache::PropagationCache) = view(prop_cache.coeffs, 1:prop_cache.active_size)
viewauxterms(prop_cache::PropagationCache) = view(prop_cache.aux_terms, 1:prop_cache.active_size)
viewauxcoeffs(prop_cache::PropagationCache) = view(prop_cache.aux_coeffs, 1:prop_cache.active_size)
viewflags(prop_cache::PropagationCache) = view(prop_cache.flags, 1:prop_cache.active_size)
viewindices(prop_cache::PropagationCache) = view(prop_cache.indices, 1:prop_cache.active_size)

term(trm::Integer) = trm
term(pstr::PauliString) = term(pstr.term)


Base.length(prop_cache::PropagationCache) = length(prop_cache.terms)
Base.isempty(prop_cache::PropagationCache) = prop_cache.active_size == 0


function topaulisum(nq, prop_cache::PropagationCache)
    TT = termtype(prop_cache)
    CT = coefftype(prop_cache)

    psum = PauliSum(nq, Dict{TT,CT}(zip(viewterms(prop_cache), viewcoeffs(prop_cache))))
    return psum
end
