## PROBLEMS:
# - We need to reuse more code from PauliPropagation.jl
# - Use truncation function that returns true or false for truncate 
# 

###
##
# Propagate a VectorPauliSum through a circuit.
# This uses the PropagationCache structure defined in PauliAlgebra/VectorPauliSum.jl
##
###

using AcceleratedKernels
const AK = AcceleratedKernels


function propagate(circuit, vecpsum::VectorPauliSum, thetas=nothing; kwargs...)
    return propagate!(circuit, deepcopy(vecpsum), thetas; kwargs...)
end


function propagate!(circuit, vecpsum::VectorPauliSum, thetas=nothing; kwargs...)
    prop_cache = PropagationCache(vecpsum)
    prop_cache = propagate!(freeze(circuit, thetas), prop_cache; kwargs...)
    vecpsum = prop_cache.vecpsum
    resize!(vecpsum, prop_cache.active_size)
    return vecpsum
end


function propagate!(circuit, prop_cache::PropagationCache; min_abs_coeff=1e-10, max_weight=Inf)
    # assume circuit contains the parameters via freezing, or is parameter-free
    @assert countparameters(circuit) == 0 "circuit requires parameters."

    for (i, gate) in enumerate(reverse(circuit))

        prop_cache = applymergetruncate!(gate, prop_cache; min_abs_coeff=min_abs_coeff, max_weight=max_weight)

    end
    return prop_cache
end

function applymergetruncate!(gate, prop_cache::PropagationCache; min_abs_coeff=1e-10, max_weight=Inf)
    prop_cache = applytoall!(gate, prop_cache)

    prop_cache = mergeterms!(prop_cache)

    prop_cache = truncate!(prop_cache; min_abs_coeff, max_weight)

    return prop_cache
end

## The apply functions

function applytoall!(gate::CliffordGate, prop_cache::PropagationCache)

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


function applytoall!(frozen_gate::FrozenGate{PauliRotation,PT}, prop_cache::PropagationCache{VT,VC,VB,VI}) where {PT,VT,VC,VB,VI}

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
    if length(prop_cache.vecpsum.terms) < n_new
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
    n_max = prop_cache.indices[prop_cache.active_size] + n

    terms_view = viewterms(prop_cache)
    coeffs = prop_cache.vecpsum.coeffs
    terms = prop_cache.vecpsum.terms
    flags = prop_cache.flags
    indices = prop_cache.indices
    @assert length(terms) >= n_max "PropagationCache terms array is not large enough to hold new terms."
    @assert length(coeffs) >= n_max "PropagationCache coeffs array is not large enough to hold new coeffs."
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

_flagstoindices!(prop_cache::PropagationCache) = _flagstoindices!(prop_cache.indices, viewflags(prop_cache))
_flagstoindices!(indices_dst, flags) = AK.accumulate!(+, indices_dst, flags; init=0)

function mergeterms!(prop_cache::PropagationCache)
    if isempty(prop_cache)
        return prop_cache
    end

    # Sort the vector to group identical keys together.
    # TODO: dispatch this in dedicated functions
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

# TODO: make this sort!(prop_cache) function instead

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
    coeffs = prop_cache.vecpsum.coeffs
    aux_terms = prop_cache.aux_vecpsum.terms
    aux_coeffs = prop_cache.aux_vecpsum.coeffs
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
    prop_cache.vecpsum, prop_cache.aux_vecpsum = prop_cache.aux_vecpsum, prop_cache.vecpsum
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
    coeffs = viewcoeffs(prop_cache)
    aux_terms = viewauxterms(prop_cache)
    aux_coeffs = viewauxcoeffs(prop_cache)
    flags = viewflags(prop_cache)
    indices = viewindices(prop_cache)
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

