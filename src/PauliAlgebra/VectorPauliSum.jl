###
##
# A file to define a Pauli sum consisting of a vector of terms and a vector of coefficients.
# Can be used for multithreaded CPU and GPU propagation.
##
###

# This type multi-threads where possible. 
using AcceleratedKernels
const AK = AcceleratedKernels


struct VectorPauliSum{TV,CV}
    nqubits::Int
    terms::TV
    coeffs::CV

    function VectorPauliSum(nqubits::Int, terms::TV, coeffs::CV) where {TV,CV}
        @assert length(terms) == length(coeffs) "Length of terms and coeffs must be the same. Got $(length(terms)) and $(length(coeffs))."
        return new{TV,CV}(nqubits, terms, coeffs)
    end
end

VectorPauliSum(nqubits::Int) = VectorPauliSum(Float64, nqubits)
VectorPauliSum(::Type{CT}, nqubits::Int) where {CT} = VectorPauliSum(nqubits, getinttype(nqubits)[], CT[])
VectorPauliSum(pstr::PauliString) = VectorPauliSum(pstr.nqubits, [pstr.term], [pstr.coeff])
VectorPauliSum(psum::PauliSum) = VectorPauliSum(psum.nqubits, collect(paulis(psum)), collect(coefficients(psum)))

paulis(vpsum::VectorPauliSum) = vpsum.terms
coefficients(vpsum::VectorPauliSum) = vpsum.coeffs
paulitype(vpsum::VectorPauliSum) = eltype(vpsum.terms)
coefftype(vpsum::VectorPauliSum) = eltype(vpsum.coeffs)

Base.similar(vpsum::VectorPauliSum) = VectorPauliSum(vpsum.nqubits, similar(vpsum.terms), similar(vpsum.coeffs))

function Base.resize!(vpsum::VectorPauliSum, n_new::Int)
    resize!(vpsum.terms, n_new)
    resize!(vpsum.coeffs, n_new)
    return vpsum
end

Base.length(vpsum::VectorPauliSum) = length(vpsum.terms)

function numcoefftype(psum::VectorPauliSum)
    if length(psum) == 0
        throw(
            "Numeric coefficient type cannot be inferred from an empty VectorPauliSum." *
            "Consider defining a `numcoefftype(psum::$(typeof(psum)))` method.")
    end
    return typeof(tonumber(first(coefficients(psum))))
end

function Base.iterate(vecpsum::VectorPauliSum)
    # 1. Create the iterator we are delegating to
    iter = zip(paulis(vecpsum), coefficients(vecpsum))

    # 2. Start its iteration
    next = iterate(iter)

    # 3. Return the first item and a new state tuple: (iterator, iterator_state)
    #    We use a ternary operator for compactness.
    return next === nothing ? nothing : (next[1], (iter, next[2]))
end

function Base.iterate(vecpsum::VectorPauliSum, state)
    # 1. Unpack the state tuple
    (iter, inner_state) = state

    # 2. Continue the delegated iteration
    next = iterate(iter, inner_state)

    # 3. Return the next item and the updated state tuple
    return next === nothing ? nothing : (next[1], (iter, next[2]))
end

function Base.show(io::IO, vecpsum::VectorPauliSum)
    n_paulis = length(vecpsum)
    if n_paulis == 0
        println(io, "Empty VectorPauliSum.")
        return
    elseif n_paulis == 1
        println(io, "VectorPauliSum with 1 term:")
    else
        println(io, "VectorPauliSum with $(n_paulis) terms:")
    end

    for i in 1:length(vecpsum)
        if i > 20
            println(io, "  ...")
            break
        end
        pauli_string = inttostring(vecpsum.terms[i], vecpsum.nqubits)
        if length(pauli_string) > 20
            pauli_string = pauli_string[1:20] * "..."
        end
        println(io, vecpsum.coeffs[i], " * $(pauli_string)")
    end
end



function Base.sort!(vpsum::VectorPauliSum; by=nothing, kwargs...)
    # instead of using sortperm, we use sort!() on an index array 
    # this is to be able to sort on any properties of the terms of coeffs 

    indices = collect(1:length(vpsum))

    # default for if "by" is not provided
    byfunc = isnothing(by) ? i -> vpsum.terms[i] : by

    AK.sort!(indices; by=byfunc, kwargs...)
    vpsum.terms .= view(vpsum.terms, indices)
    vpsum.coeffs .= view(vpsum.coeffs, indices)
    return vpsum
end



###
##
# A PropagationCache carries two VectorPauliSums and flags and indices for propagation.
# It is used inside propagate() to avoid reallocations.
##
###

mutable struct PropagationCache{VT,VC,VB,VI}
    vecpsum::VectorPauliSum{VT,VC}
    aux_vecpsum::VectorPauliSum{VT,VC}
    flags::VB
    indices::VI

    # we will over-allocate the arrays and keep track of the non-empty size
    active_size::Int
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
        pauli_string = inttostring(prop_cache.vecpsum.terms[i], prop_cache.vecpsum.nqubits)
        if length(pauli_string) > 20
            pauli_string = pauli_string[1:20] * "..."
        end
        println(io, prop_cache.coeffs[i], " * $(pauli_string)")
    end
end


function Base.resize!(prop_cache::PropagationCache, n_new::Int)
    resize!(prop_cache.vecpsum, n_new)
    resize!(prop_cache.aux_vecpsum, n_new)
    resize!(prop_cache.flags, n_new)
    resize!(prop_cache.indices, n_new)
    return prop_cache
end

PauliPropagation.paulitype(prop_cache::PropagationCache{VT,VC,VB,VI}) where {VT,VC,VB,VI} = eltype(VT)
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
