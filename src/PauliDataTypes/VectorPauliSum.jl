###
##
# A file to define a Pauli sum consisting of a vector of terms and a vector of coefficients.
# Can be used for multithreaded CPU and GPU propagation.
##
###

# This type multi-threads where possible. 
using AcceleratedKernels
const AK = AcceleratedKernels


struct VectorPauliSum{TV,CV} <: AbstractPauliSum
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


"""
    nqubits(vpsum::VectorPauliSum)

Get the number of qubits that the `VectorPauliSum` is defined on.
"""
nqubits(vpsum::VectorPauliSum) = vpsum.nqubits

"""
    paulis(vpsum::VectorPauliSum)

Get the vector of Pauli strings from a `VectorPauliSum`.
"""
paulis(vpsum::VectorPauliSum) = vpsum.terms

"""
    coefficients(vpsum::VectorPauliSum)
    
Get the vector of coefficients from a `VectorPauliSum`.
"""
coefficients(vpsum::VectorPauliSum) = vpsum.coeffs

"""
    paulitype(vpsum::VectorPauliSum)

Get the Pauli integer type of a `VectorPauliSum`.
"""
paulitype(vpsum::VectorPauliSum) = eltype(vpsum.terms)

"""
    coefftype(vpsum::VectorPauliSum)

Get the coefficient type of a `VectorPauliSum`.
"""
coefftype(vpsum::VectorPauliSum) = eltype(vpsum.coeffs)


"""
    getcoeff(vpsum::VectorPauliSum, pstr::Integer)

Get the coefficient of an integer Pauli strings `pstr` in a `PauliSum`. 
Defaults to 0 if the Pauli string is not in the `PauliSum`.
"""
function getcoeff(vpsum::VectorPauliSum, pstr::Integer)
    return sum(vpsum.coeffs[i] for i in eachindex(vpsum.terms) if vpsum.terms[i] == pstr; init=zero(coefftype(vpsum)))
end


"""
    topaulistrings(vpsum::VectorPauliSum)

Returns the Pauli strings in a `PauliSum` and their coefficients as a list of `PauliString`.
"""
topaulistrings(vpsum::VectorPauliSum) = [PauliString(vpsum.nqubits, pauli, coeff) for (pauli, coeff) in zip(vpsum.terms, vpsum.coeffs)]


Base.similar(vpsum::VectorPauliSum) = VectorPauliSum(vpsum.nqubits, similar(vpsum.terms), similar(vpsum.coeffs))

function Base.resize!(vpsum::VectorPauliSum, n_new::Int)
    resize!(vpsum.terms, n_new)
    resize!(vpsum.coeffs, n_new)
    return vpsum
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

function PropagationCache(vpsum::PauliSum)
    return PropagationCache(VectorPauliSum(vpsum))
end

# convert back
function VectorPauliSum(prop_cache::PropagationCache)
    vecpsum = deepcopy(prop_cache.vecpsum)
    resize!(vecpsum, prop_cache.active_size)
    return vecpsum
end

function PauliSum(prop_cache::PropagationCache)
    mergeterms!(prop_cache)
    return PauliSum(prop_cache.vecpsum.nqubits, Dict(zip(viewterms(prop_cache), viewcoeffs(prop_cache))))
end

"""
    nqubits(prop_cache::PropagationCache)

Get the number of qubits that the `PropagationCache` is defined on.
"""
nqubits(prop_cache::PropagationCache) = prop_cache.vecpsum.nqubits


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
        println(io, prop_cache.vecpsum.coeffs[i], " * $(pauli_string)")
    end
end


function Base.resize!(prop_cache::PropagationCache, n_new::Int)
    resize!(prop_cache.vecpsum, n_new)
    resize!(prop_cache.aux_vecpsum, n_new)
    resize!(prop_cache.flags, n_new)
    resize!(prop_cache.indices, n_new)
    return prop_cache
end


paulis(prop_cache) = viewterms(prop_cache)
coefficients(prop_cache) = viewcoeffs(prop_cache)
paulitype(prop_cache::PropagationCache{VT,VC,VB,VI}) where {VT,VC,VB,VI} = eltype(VT)
coefftype(prop_cache::PropagationCache{VT,VC,VB,VI}) where {VT,VC,VB,VI} = eltype(VC)

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

