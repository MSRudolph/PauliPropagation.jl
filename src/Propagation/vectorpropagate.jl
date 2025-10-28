using Base.Threads
using AcceleratedKernels
const AK = AcceleratedKernels
using PauliPropagation


mutable struct PropagationCache{VT,VC,VB,VI}
    terms::VT
    coeffs::VC
    aux_terms::VT
    aux_coeffs::VC
    flags::VB
    indices::VI

    # we will over-allocate the arrays and keep track of the non-empty size
    active_size::Int
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

function PropagationCache(terms, coeffs)
    aux_terms = similar(terms)
    aux_coeffs = similar(coeffs)
    flags = similar(terms, Bool)
    indices = similar(terms, Int)
    return PropagationCache(terms, coeffs, aux_terms, aux_coeffs, flags, indices, length(terms))
end

function PropagationCache(pstrings::AbstractVector{PauliString{TT,CT}}) where {TT,CT}
    terms = similar(pstrings, TT)
    coeffs = similar(pstrings, CT)
    for (i, pstr) in enumerate(pstrings)
        terms[i] = pstr.term
        coeffs[i] = pstr.coeff
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
        println(io, prop_cache.coeffs[i], " * $(prop_cache.terms[i])")
    end
end

termtype(prop_cache::PropagationCache{VT,VC,VB,VI}) where {VT,VC,VB,VI} = eltype(VT)
coefftype(prop_cache::PropagationCache{VT,VC,VB,VI}) where {VT,VC,VB,VI} = eltype(VC)

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


function vectorpropagate(circuit, strings::AbstractVector{<:PauliString}, args...; kwargs...)
    return vectorpropagate(circuit, PropagationCache(strings), args...; kwargs...).terms
end

function vectorpropagate(circuit, prop_cache::PropagationCache, thetas; kwargs...)
    return vectorpropagate(freeze(circuit, thetas), prop_cache; kwargs...)
end

function vectorpropagate(circuit, prop_cache::PropagationCache; min_abs_coeff=1e-10, max_weight=Inf, kwargs...)
    # assume circuit contains the parameters via freezing, or is parameter-free
    @assert countparameters(circuit) == 0 "'circuit' must be parameter-free. Consider using 'freeze()'."

    for (i, gate) in enumerate(reverse(circuit))

        prop_cache = applygate!(gate, prop_cache)

        prop_cache = mergeterms!(prop_cache)

        prop_cache = truncate!(prop_cache; min_abs_coeff, max_weight)

    end
    return prop_cache
end

## The apply functions

function applygate!(gate::CliffordGate, prop_cache::PropagationCache)

    # everything is done in place
    terms_view = viewterms(prop_cache)
    coeffs_view = viewcoeffs(prop_cache)
    @assert length(terms_view) == length(coeffs_view)
    AK.foreachindex(terms_view) do ii
        term = terms_view[ii]
        coeff = coeffs_view[ii]

        term, coeff = apply(gate, term, coeff)

        # inbounds is safe here because we assert equal lengths
        terms_view[ii] = term
        coeffs_view[ii] = coeff
    end

    return prop_cache
end


function _paulirotationproduct(gate_mask::TT, pstr::TT) where TT
    new_pstr = PauliPropagation._bitpaulimultiply(gate_mask, pstr)

    # this counts the exponent of the imaginary unit in the new Pauli string
    im_count = PauliPropagation._calculatesignexponent(gate_mask, pstr)

    # now, instead of computing im^im_count followed by another im factor from the gate rules,
    # we do this in one step via a cheeky trick:
    sign = (im_count & 2) - 1
    # this is equivalent to sign = real( im * im^im_count)

    return new_pstr, sign
end

function applygate!(frozen_gate::FrozenGate{PauliRotation,PT}, prop_cache::PropagationCache{VT,VC,VB,VI}) where {PT,VT,VC,VB,VI}

    TT = eltype(VT)

    n_old = prop_cache.active_size

    # unpack the frozen PauliRotation
    gate, theta = frozen_gate.gate, frozen_gate.parameter

    # this allows for faster operations
    masked_gate = PauliPropagation._tomaskedpaulirotation(gate, TT)

    # get the mask out because because the gate cannot be in the function when using GPU
    gate_mask = masked_gate.generator_mask

    # this needs to be in a separate function because variable names cannot be duplicated (WOW)
    _flaganticommuting!(prop_cache, gate_mask)

    n_noncommutes = prop_cache.indices[prop_cache.active_size]

    # slit off into the same array
    n_new = n_old + n_noncommutes

    # potential resize factor
    resize_factor = 2
    if length(prop_cache.terms) < n_new
        resize!(prop_cache, n_new * resize_factor)
    end

    _applypaulirotation!(prop_cache, gate_mask, theta)
    prop_cache.active_size = n_new

    return prop_cache
end

function _applypaulirotation!(prop_cache::PropagationCache, gate_mask::TT, theta) where {TT}

    # pre-compute the sine and cosine values because the are used for every Pauli string that does not commute with the gate
    cos_val = cos(theta)
    sin_val = sin(theta)

    n = prop_cache.active_size

    terms_view = viewterms(prop_cache)
    coeffs = prop_cache.coeffs
    terms = prop_cache.terms
    flags = prop_cache.flags
    indices = prop_cache.indices
    AK.foreachindex(terms_view) do ii
        # here it anticommutes
        if flags[ii]
            term = terms[ii]
            coeff = coeffs[ii]

            coeff1 = coeff * cos_val
            new_term, sign = _paulirotationproduct(gate_mask, term)
            coeff2 = coeff * sin_val * sign

            coeffs[ii] = coeff1

            terms[n+indices[ii]] = new_term
            coeffs[n+indices[ii]] = coeff2
        end
    end

    return
end

function _flaganticommuting!(prop_cache::PropagationCache, gate_mask::TT) where {TT}
    terms_view = viewterms(prop_cache)
    flags = prop_cache.flags
    AK.foreachindex(terms_view) do ii
        flags[ii] = !commutes(terms_view[ii], gate_mask)
    end

    _flagstoindices!(prop_cache)
    return
end

_flagstoindices!(prop_cache::PropagationCache) = _flagstoindices!(prop_cache.indices, viewflags(prop_cache))
_flagstoindices!(indices_dst, flags) = AK.accumulate!(+, indices_dst, flags; init=0)

function mergeterms!(prop_cache::PropagationCache)
    if isempty(prop_cache)
        return prop_cache
    end

    # Sort the vector to group identical keys together.
    AK.sortperm!(viewindices(prop_cache), viewterms(prop_cache), by=term)

    # shuffle the terms and coeffs according to the sorted indices into their aux arrays
    _permuteviaindices!(prop_cache)
    swapterms!(prop_cache)

    # Find the start of each group in parallel.
    _flagunique!(prop_cache)

    # get the indices for where to add the coefficients
    _flagstoindices!(prop_cache)
    # TODO: finalindex() function to work on GPU and CPU without sharing memory
    n_unique_terms = prop_cache.indices[prop_cache.active_size]

    # early stop if all are unique 
    if n_unique_terms == prop_cache.active_size
        return prop_cache
    end

    # the reason we don't do atomic add instead is because the entry might not start off as 0
    _deduplicate!(prop_cache)

    prop_cache.active_size = n_unique_terms

    return prop_cache
end

function _permuteviaindices!(prop_cache::PropagationCache)
    indices_view = viewindices(prop_cache)
    term_view = viewterms(prop_cache)
    coeffs_view = viewcoeffs(prop_cache)
    aux_terms_view = viewauxterms(prop_cache)
    aux_coeffs_view = viewauxcoeffs(prop_cache)
    AK.foreachindex(indices_view) do ii
        sorted_idx = indices_view[ii]
        aux_terms_view[ii] = term_view[sorted_idx]
        aux_coeffs_view[ii] = coeffs_view[sorted_idx]
    end
    return
end

function _flagunique!(prop_cache::PropagationCache)
    term_view = viewterms(prop_cache)
    flags_view = viewflags(prop_cache)
    AK.foreachindex(term_view) do ii
        if ii == 1
            flags_view[ii] = true
        else
            flags_view[ii] = term_view[ii] != term_view[ii-1]
        end
    end
    return
end

function _deduplicate!(prop_cache::PropagationCache)

    term_view = viewterms(prop_cache)
    coeffs = prop_cache.coeffs
    aux_terms = prop_cache.aux_terms
    aux_coeffs = prop_cache.aux_coeffs
    flags = prop_cache.flags
    indices = prop_cache.indices
    active_size = prop_cache.active_size
    AK.foreachindex(term_view) do ii
        # if this is the start of a new group
        if flags[ii]
            # end index is the before the next flag or the end of the array
            end_idx = ii
            while end_idx < active_size && !flags[end_idx+1]
                end_idx += 1
            end

            # Sum the values in the range.
            CT = typeof(coeffs[ii])
            merged_coeff = zero(CT)
            for jj in ii:end_idx
                merged_coeff += coeffs[jj]
            end

            aux_terms[indices[ii]] = term_view[ii]
            aux_coeffs[indices[ii]] = merged_coeff
        end
    end

    # swap terms and aux_terms
    swapterms!(prop_cache)

    return
end

function swapterms!(prop_cache::PropagationCache)
    prop_cache.terms, prop_cache.aux_terms = prop_cache.aux_terms, prop_cache.terms
    prop_cache.coeffs, prop_cache.aux_coeffs = prop_cache.aux_coeffs, prop_cache.coeffs
    return prop_cache
end


function truncate!(prop_cache::PropagationCache{TT,CT}; min_abs_coeff::Real, max_weight::Real=Inf) where {TT,CT}

    if isempty(prop_cache)
        return prop_cache
    end

    # TODO: can we simplify this entire function via a parallel filter!() function?

    # flag the indices that we keep
    _flagtokeep!(prop_cache, min_abs_coeff, max_weight)

    # get the new indices after deletion
    _flagstoindices!(prop_cache)
    n_kept = prop_cache.indices[prop_cache.active_size]

    _moveflagged!(prop_cache)

    swapterms!(prop_cache)
    prop_cache.active_size = n_kept

    return prop_cache
end

function _moveflagged!(prop_cache::PropagationCache{TT,CT}) where {TT,CT}
    terms_view = viewterms(prop_cache)
    coeffs = prop_cache.coeffs
    aux_terms = prop_cache.aux_terms
    aux_coeffs = prop_cache.aux_coeffs
    flags = prop_cache.flags
    indices = prop_cache.indices
    AK.foreachindex(terms_view) do ii
        if flags[ii]
            aux_terms[indices[ii]] = terms_view[ii]
            aux_coeffs[indices[ii]] = coeffs[ii]
        end
    end
    return
end

function _flagtokeep!(prop_cache::PropagationCache{TT,CT}, min_abs_coeff, max_weight) where {TT,CT}
    terms_view = viewterms(prop_cache)
    coeffs_view = viewcoeffs(prop_cache)
    @assert length(terms_view) == length(coeffs_view)
    flags = prop_cache.flags
    AK.foreachindex(terms_view) do ii
        flags[ii] = (abs(tonumber(coeffs_view[ii])) >= min_abs_coeff) && countweight(terms_view[ii]) <= max_weight
    end
    return
end

# function chopterms!(prop_cache::PropagationCache, max_terms)
#     if length(prop_cache.terms) <= max_terms
#         return prop_cache
#     end

#     # sort by coefficient magnitude
#     AcceleratedKernels.sort!(
#         prop_cache.terms,
#         by=p -> abs(tonumber(p.coeff));
#         rev=true
#     )

#     # keep only the largest 'max_terms' terms
#     resize!(prop_cache.terms, max_terms)
#     resize!(prop_cache.flags, max_terms)
#     resize!(prop_cache.indices, max_terms)
#     return prop_cache
# end

function zeroinactive!(prop_cache::PropagationCache)
    TT = termtype(prop_cache)
    CT = coefftype(prop_cache)

    @assert length(prop_cache.terms) == length(prop_cache.coeffs)

    terms = prop_cache.terms
    coeffs = prop_cache.coeffs
    active_size = prop_cache.active_size
    AK.foreachindex(terms) do ii
        if ii > active_size
            terms[ii] = zero(TT)
            coeffs[ii] = zero(CT)
        end
    end
    return prop_cache
end

function topaulisum(nq, prop_cache::PropagationCache)
    TT = termtype(prop_cache)
    CT = coefftype(prop_cache)

    psum = PauliSum(nq, Dict{TT,CT}(zip(viewterms(prop_cache), viewcoeffs(prop_cache))))
    return psum
end


function PauliPropagation.overlapwithzero(strings::Vector)
    total = zero(numcoefftype(strings[1]))
    for pstr in strings

        total += tonumber(pstr.coeff) * !containsXorY(pstr.term)
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