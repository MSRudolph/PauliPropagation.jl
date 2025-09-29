using Base.Threads
using AcceleratedKernels
using PauliPropagation

struct WeightedString{TT,CT}
    term::TT
    coeff::CT
end


# TODO: an is_sorted flag?
struct PropagationCache{TT,CT}
    terms::Vector{WeightedString{TT,CT}}
    flags::Vector{Bool}
    indices::Vector{Int}

    function PropagationCache(terms::Vector{WeightedString{TT,CT}}, flags, indices) where {TT,CT}
        @assert length(terms) == length(flags) == length(indices)
        new{TT,CT}(terms, flags, indices)
    end
end

function PropagationCache(strings::Vector)
    n = length(strings)
    flags = Vector{Bool}(undef, n)
    indices = Vector{Int}(undef, n)
    return PropagationCache(strings, flags, indices)
end


function vectorpropagate(circuit, strings::Vector; kwargs...)
    return vectorpropagate(circuit, PropagationCache(strings); kwargs...).terms
end

function vectorpropagate(circuit, prop_cache::PropagationCache; min_abs_coeff=1e-10, kwargs...)
    # assume circuit contains the parameters via freezing, or is parameter-free
    @assert countparameters(circuit) == 0 "'circuit' must be parameter-free. Consider using 'freeze()'."

    n_gates = length(circuit)
    for (i, gate) in enumerate(reverse(circuit))

        prop_cache = applygate!(gate, prop_cache)

        prop_cache = mergeterms!(prop_cache)
        prop_cache = truncate!(prop_cache, min_abs_coeff)

        ## Thoughts:
        # - Can we merge and truncate in one?
        # - Is coeff truncation strictly worse than max cache size?
        # - Can we write a sorting function between the two sorted arrays?


    end
    return prop_cache
end

## The apply functions

function applygate!(gate, pstrings::Vector)
    prop_cache = PropagationCache(pstrings)
    prop_cache = applygate!(gate, prop_cache)
    return prop_cache.terms

end

function applygate!(gate::CliffordGate, prop_cache::PropagationCache)

    # everything is done in place
    @threads for ii in eachindex(prop_cache.terms)
        wstring = prop_cache.terms[ii]

        term, coeff = apply(gate, wstring.term, wstring.coeff)

        prop_cache.terms[ii] = WeightedString(term, coeff)
    end

    return prop_cache
end

function applygate!(frozen_gate::FrozenGate{PauliRotation,PT}, prop_cache::PropagationCache{TT,CT}) where {CT,TT,PT}

    n_oldstrings = length(prop_cache.terms)

    # unpack the frozen PauliRotation
    gate, theta = frozen_gate.gate, frozen_gate.parameter

    # this allows for faster operations
    masked_gate = PauliPropagation._tomaskedpaulirotation(gate, TT)

    # pre-compute the sine and cosine values because the are used for every Pauli string that does not commute with the gate
    cos_val = cos(theta)
    sin_val = sin(theta)

    @threads for ii in eachindex(prop_cache.terms)
        wstring = prop_cache.terms[ii]
        prop_cache.flags[ii] = !commutes(masked_gate, wstring.term)
    end

    # index_vector = cumsum(anticommutes_vector)
    AcceleratedKernels.accumulate!(+, prop_cache.indices, prop_cache.flags; init=0)
    n_noncommutes = prop_cache.indices[end]

    # pre-allocate the output array
    new_pauli_strings = Vector{WeightedString{TT,CT}}(undef, n_oldstrings + n_noncommutes)
    # new_prop_cache = PropagationCache(new_pauli_strings)

    @threads for ii in eachindex(prop_cache.terms)
        wstring = prop_cache.terms[ii]

        if prop_cache.flags[ii]
            term = wstring.term
            coeff = wstring.coeff

            coeff1 = coeff * cos_val
            new_term, sign = PauliPropagation.getnewpaulistring(masked_gate, term)
            coeff2 = coeff * sin_val * sign

            new_pauli_strings[ii] = WeightedString(term, coeff1)
            new_pauli_strings[n_oldstrings+prop_cache.indices[ii]] = WeightedString(new_term, coeff2)

        else
            new_pauli_strings[ii] = wstring
        end
    end

    # re-use the flags and indices arrays
    resize!(prop_cache.flags, length(new_pauli_strings))
    resize!(prop_cache.indices, length(new_pauli_strings))
    return PropagationCache(new_pauli_strings, prop_cache.flags, prop_cache.indices)
end

term(wstring::WeightedString) = wstring.term

function mergeterms!(prop_cache::PropagationCache{TT,CT}) where {TT,CT}
    if isempty(prop_cache.terms)
        return prop_cache
    end

    n_terms = length(prop_cache.terms)

    # Sort the vector to group identical keys together.
    AcceleratedKernels.sort!(prop_cache.terms, by=term)

    # Find the start of each group in parallel.
    prop_cache.flags[1] = true # The first element always starts a group.
    @threads for i in 2:n_terms
        # A new group starts if the term is different from the previous one.
        prop_cache.flags[i] = prop_cache.terms[i].term != prop_cache.terms[i-1].term
    end

    # get the indices for where to add the coefficients
    AcceleratedKernels.accumulate!(+, prop_cache.indices, prop_cache.flags; init=0)
    n_unique_terms = prop_cache.indices[end]

    # create the new cache
    merged_prop_cache = PropagationCache(typeof(prop_cache.terms)(undef, n_unique_terms))

    # populate the start_indices vector
    @threads for i in eachindex(prop_cache.terms)
        if prop_cache.flags[i]
            # `prop_cache.indices[i]` gives the 1-based index for this new group.
            # We place the actual array index `i` into that slot.
            merged_prop_cache.indices[prop_cache.indices[i]] = i
        end
    end

    # sum over each group
    @threads for i in 1:n_unique_terms
        start_idx = merged_prop_cache.indices[i]
        # The group ends right before the next group starts, or at the end of the array.
        end_idx = (i < n_unique_terms) ? merged_prop_cache.indices[i+1] - 1 : n_terms

        # The term is the same for all elements in this range.
        group_term = prop_cache.terms[start_idx].term

        # Sum the values in the range.
        total_sum = zero(CT)

        for j in start_idx:end_idx
            total_sum += prop_cache.terms[j].coeff
        end

        merged_prop_cache.terms[i] = WeightedString(group_term, total_sum)
    end

    return merged_prop_cache
end


function truncate!(prop_cache::PropagationCache{TT,CT}, min_abs_coeff::Real) where {TT,CT}

    # flag the qubits that would not be deleted
    @threads for ii in eachindex(prop_cache.terms)
        pstr = prop_cache.terms[ii]
        prop_cache.flags[ii] = abs(tonumber(pstr.coeff)) >= min_abs_coeff
    end

    # get the new indices after deletion
    AcceleratedKernels.accumulate!(+, prop_cache.indices, prop_cache.flags; init=0)
    n_kept = prop_cache.indices[end]

    # pre-allocate the output array
    filtered_terms = Vector{WeightedString{TT,CT}}(undef, n_kept)
    @threads for ii in eachindex(prop_cache.terms)
        if prop_cache.flags[ii]
            filtered_terms[prop_cache.indices[ii]] = prop_cache.terms[ii]
        end
    end
    # resize old arrays 
    flags = resize!(prop_cache.flags, n_kept)
    indices = resize!(prop_cache.indices, n_kept)
    return PropagationCache(filtered_terms, flags, indices)
end

function chopterms!(prop_cache::PropagationCache, max_terms)
    if length(prop_cache.terms) <= max_terms
        return prop_cache
    end

    # sort by coefficient magnitude
    AcceleratedKernels.sort!(
        prop_cache.terms,
        by=p -> abs(tonumber(p.coeff));
        rev=true
    )

    # keep only the largest 'max_terms' terms
    resize!(prop_cache.terms, max_terms)
    resize!(prop_cache.flags, max_terms)
    resize!(prop_cache.indices, max_terms)
    return prop_cache
end



function PauliPropagation.overlapwithzero(strings::Vector)
    total = zero(numcoefftype(strings[1]))
    for pstr in strings

        total += tonumber(pstr.coeff) * !containsXorY(pstr.term)
    end
    return total
end


function PauliPropagation.overlapwithcomputational(strings::Vector{WeightedString{TT,CT}}, onebitinds) where {TT,CT}
    val = zero(CT)
    for wstring in strings
        val += tonumber(wstring.coeff) * PauliPropagation._calcsignwithones(wstring.term, onebitinds)
    end
    return val
end