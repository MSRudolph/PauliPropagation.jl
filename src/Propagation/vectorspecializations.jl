## PROBLEMS:
# - We need to reuse more code from PauliPropagation.jl
# - Use truncation function that returns true or false for truncate 
# - Implement symmetry merging 
# - Implement state overlaps
# - Add docstrings
# - Add tests
# 

###
##
# Propagate a VectorPauliSum through a circuit.
# This uses the VectorPauliPropagationCache structure defined in PauliAlgebra/VectorPauliSum.jl
##
###


# function propagate(circuit, vecpsum::VectorPauliSum, thetas=nothing; kwargs...)
#     return propagate!(circuit, deepcopy(vecpsum), thetas; kwargs...)
# end


function PropagationBase.propagate!(circuit, vecpsum::VectorPauliSum, args...; kwargs...)
    prop_cache = VectorPauliPropagationCache(vecpsum)
    prop_cache = propagate!(circuit, prop_cache, args...; kwargs...)
    vecpsum = mainsum(prop_cache)
    resize!(vecpsum, activesize(prop_cache))
    return vecpsum
end


# function propagate!(circuit, prop_cache::VectorPauliPropagationCache; min_abs_coeff=1e-10, max_weight=Inf)
#     # assume circuit contains the parameters via freezing, or is parameter-free
#     @assert countparameters(circuit) == 0 "circuit requires parameters."

#     for (i, gate) in enumerate(reverse(circuit))

#         prop_cache = applymergetruncate!(gate, prop_cache; min_abs_coeff=min_abs_coeff, max_weight=max_weight)

#     end
#     return prop_cache
# end

# function applymergetruncate!(gate, prop_cache::VectorPauliPropagationCache; min_abs_coeff=1e-10, max_weight=Inf)
#     prop_cache = applytoall!(gate, prop_cache)

#     prop_cache = mergeterms!(prop_cache)

#     prop_cache = truncate!(prop_cache; min_abs_coeff, max_weight)

#     return prop_cache
# end

## The apply functions
# TODO: overload applymergetruncate!() and don't merge
function PropagationBase.applytoall!(gate::CliffordGate, prop_cache::VectorPauliPropagationCache)
    # TODO: This needs to be reworked for GPU support

    # everything is done in place
    terms_view = activeterms(prop_cache)
    coeffs_view = activecoeffs(prop_cache)
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


function PropagationBase.applytoall!(gate::PauliRotation, prop_cache::VectorPauliPropagationCache, theta)

    # TODO: design this function in a way that it can be the default for branching gates. 
    # Think of U3 or amplitude damping 

    if prop_cache.active_size == 0
        return prop_cache
    end

    n_old = prop_cache.active_size

    # get the mask out because because the gate cannot be in the function when using GPU
    gate_mask = symboltoint(nqubits(prop_cache), gate.symbols, gate.qinds)

    # this needs to be in a separate function because variable names cannot be duplicated (WOW)
    # _flaganticommuting!(prop_cache, gate_mask)
    flagterms!(trm -> !commutes(trm, gate_mask), prop_cache)

    # this runs a cumsum over the flags to get the indices
    flagstoindices!(prop_cache)

    # the final index is the number of new terms
    n_noncommutes = lastactiveindex(prop_cache)

    # slit off into the same array
    n_new = n_old + n_noncommutes

    # potential resize factor
    resize_factor = 2
    if capacity(prop_cache) < n_new
        resize!(prop_cache, n_new * resize_factor)
    end

    # does the branching logic
    _applypaulirotation!(prop_cache, gate_mask, theta)

    # we now have n_new possibly douplicate Pauli strings in the array
    setactivesize!(prop_cache, n_new)

    return prop_cache
end

function _applypaulirotation!(prop_cache::VectorPauliPropagationCache, gate_mask::TT, theta) where {TT}

    # pre-compute the sine and cosine values because the are used for every Pauli string that does not commute with the gate
    cos_val = cos(theta)
    sin_val = sin(theta)

    n = activesize(prop_cache)
    n_max = n + lastactiveindex(prop_cache)

    active_terms = activeterms(prop_cache)

    # full-length terms so we can write new terms at the end
    terms = paulis(mainsum(prop_cache))
    coeffs = coefficients(mainsum(prop_cache))
    @assert length(terms) >= n_max "VectorPauliPropagationCache terms array is not large enough to hold new terms."
    @assert length(coeffs) >= n_max "VectorPauliPropagationCache coeffs array is not large enough to hold new coeffs."

    flags = activeflags(prop_cache)
    indices = activeindices(prop_cache)

    # TODO: modularize this into something like "two-branching pattern"
    AK.foreachindex(active_terms) do ii
        # here it anticommutes
        if flags[ii]
            term = terms[ii]
            coeff = coeffs[ii]

            coeff1 = coeff * cos_val
            new_term, sign = paulirotationproduct(gate_mask, term)
            coeff2 = coeff * sin_val * sign

            coeffs[ii] = coeff1

            terms[n+indices[ii]] = new_term
            coeffs[n+indices[ii]] = coeff2
        end
    end

    return
end


# function mergeterms!(prop_cache::VectorPauliPropagationCache)
#     if isempty(prop_cache)
#         return prop_cache
#     end

#     # Sort the vector to group identical keys together.
#     # TODO: dispatch this in dedicated functions
#     AK.sortperm!(viewindices(prop_cache), viewterms(prop_cache), by=term)

#     # shuffle the terms and coeffs according to the sorted indices into their aux arrays
#     permuteviaindices!(prop_cache)

#     # Find the start of each group in parallel.
#     _flagunique!(prop_cache)

#     # get the indices for where to add the coefficients
#     flagstoindices!(prop_cache)
#     # TODO: finalindex() function to work on GPU and CPU without sharing memory
#     n_unique_terms = prop_cache.indices[prop_cache.active_size]

#     # early stop if all are unique 
#     if n_unique_terms == prop_cache.active_size
#         return prop_cache
#     end

#     # the reason we don't do atomic add instead is because the entry might not start off as 0
#     _deduplicate!(prop_cache)

#     prop_cache.active_size = n_unique_terms

#     return prop_cache
# end




# function _flagunique!(prop_cache::VectorPauliPropagationCache)
#     term_view = activeterms(prop_cache)
#     flags_view = activeflags(prop_cache)
#     AK.foreachindex(term_view) do ii
#         if ii == 1
#             flags_view[ii] = true
#         else
#             flags_view[ii] = term_view[ii] != term_view[ii-1]
#         end
#     end
#     return prop_cache
# end

# function _deduplicate!(prop_cache::VectorPauliPropagationCache)

#     term_view = activeterms(prop_cache)
#     coeffs = prop_cache.psum.coeffs
#     aux_terms = prop_cache.aux_psum.terms
#     aux_coeffs = prop_cache.aux_psum.coeffs
#     flags = prop_cache.flags
#     indices = prop_cache.indices
#     active_size = prop_cache.active_size
#     AK.foreachindex(term_view) do ii
#         # if this is the start of a new group
#         if flags[ii]
#             # end index is the before the next flag or the end of the array
#             end_idx = ii
#             while end_idx < active_size && !flags[end_idx+1]
#                 end_idx += 1
#             end

#             # Sum the values in the range.
#             CT = typeof(coeffs[ii])
#             merged_coeff = zero(CT)
#             for jj in ii:end_idx
#                 merged_coeff += coeffs[jj]
#             end

#             aux_terms[indices[ii]] = term_view[ii]
#             aux_coeffs[indices[ii]] = merged_coeff
#         end
#     end

#     # swap terms and aux_terms
#     swapsums!(prop_cache)

#     return prop_cache
# end


# function truncate!(prop_cache::VectorPauliPropagationCache{TT,CT}; min_abs_coeff::Real, max_weight::Real=Inf) where {TT,CT}
#     # TODO: in Base this should take min_abs_coeff and a custom truncation function. 
#     # Specialized basis libraries overload truncate!() for their cache type,
#     # but just use the function to load the custom truncation function to pass to Base.
#     if isempty(prop_cache)
#         return prop_cache
#     end

#     # TODO: can we simplify this entire function via a parallel filter!() function?

#     # flag the indices that we keep
#     _flagtokeep!(prop_cache, min_abs_coeff, max_weight)

#     # get the new indices after deletion
#     flagstoindices!(prop_cache)

#     _moveflagged!(prop_cache)

#     swapsums!(prop_cache)

#     n_kept = prop_cache.indices[prop_cache.active_size]
#     prop_cache.active_size = n_kept

#     return prop_cache
# end

# function filterviaflags!(prop_cache::VectorPauliPropagationCache{TT,CT}) where {TT,CT}
#     terms_view = viewterms(prop_cache)
#     coeffs = viewcoeffs(prop_cache)
#     aux_terms = viewauxterms(prop_cache)
#     aux_coeffs = viewauxcoeffs(prop_cache)
#     flags = viewflags(prop_cache)
#     indices = viewindices(prop_cache)

#     filterviaflags!(aux_terms, aux_coeffs, terms_view, coeffs, flags, indices)

#     swapsums!(prop_cache)

#     n_new = indices[prop_cache.active_size]
#     prop_cache.active_size = n_new

#     return prop_cache
# end

# function _moveflagged!(prop_cache::VectorPauliPropagationCache{TT,CT}) where {TT,CT}
#     terms_view = activeterms(prop_cache)
#     coeffs = activecoeffs(prop_cache)
#     aux_terms = activeauxterms(prop_cache)
#     aux_coeffs = activeauxcoeffs(prop_cache)
#     flags = activeflags(prop_cache)
#     indices = activeindices(prop_cache)
#     AK.foreachindex(terms_view) do ii
#         if flags[ii]
#             aux_terms[indices[ii]] = terms_view[ii]
#             aux_coeffs[indices[ii]] = coeffs[ii]
#         end
#     end
#     return
# end

# function _flagtokeep!(prop_cache::VectorPauliPropagationCache{TT,CT}, min_abs_coeff, max_weight) where {TT,CT}
#     terms_view = activeterms(prop_cache)
#     coeffs_view = activecoeffs(prop_cache)
#     @assert length(terms_view) == length(coeffs_view)
#     flags = prop_cache.flags
#     AK.foreachindex(terms_view) do ii
#         flags[ii] = (abs(tonumber(coeffs_view[ii])) >= min_abs_coeff) && countweight(terms_view[ii]) <= max_weight
#     end
#     return
# end
