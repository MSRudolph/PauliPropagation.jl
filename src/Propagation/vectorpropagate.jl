using Base.Threads
using AcceleratedKernels
const AK = AcceleratedKernels
using PauliPropagation

struct WeightedString{TT,CT}
    term::TT
    coeff::CT
end

termtype(wstring::WeightedString{TT,CT}) where {TT,CT} = TT
coefftype(wstring::WeightedString{TT,CT}) where {TT,CT} = CT

Base.:*(wstring::WeightedString, c) = WeightedString(wstring.term, wstring.coeff * c)

# types are:
# vector of terms 
# vector of flags (bool)
# vector of indices (int)
mutable struct PropagationCache{VT,VB,VI}
    terms::VT
    aux_terms::VT
    flags::VB
    indices::VI

    # we will over-allocate the arrays and keep track of the non-empty size
    active_size::Int
end

function Base.resize!(prop_cache::PropagationCache, n_new::Int)
    resize!(prop_cache.terms, n_new)
    resize!(prop_cache.aux_terms, n_new)
    resize!(prop_cache.flags, n_new)
    resize!(prop_cache.indices, n_new)
    return prop_cache
end

function PropagationCache(terms)
    aux_terms = similar(terms)
    flags = similar(terms, Bool)
    indices = similar(terms, Int)
    return PropagationCache(terms, aux_terms, flags, indices, length(terms))
end

termtype(prop_cache::PropagationCache{VT,VB,VI}) where {VT,VB,VI} = eltype(VT).parameters[1]
coefftype(prop_cache::PropagationCache{VT,VB,VI}) where {VT,VB,VI} = eltype(VT).parameters[2]

viewterms(prop_cache::PropagationCache) = view(prop_cache.terms, 1:prop_cache.active_size)
viewauxterms(prop_cache::PropagationCache) = view(prop_cache.aux_terms, 1:prop_cache.active_size)
viewflags(prop_cache::PropagationCache) = view(prop_cache.flags, 1:prop_cache.active_size)
viewindices(prop_cache::PropagationCache) = view(prop_cache.indices, 1:prop_cache.active_size)

Base.length(prop_cache::PropagationCache) = prop_cache.active_size
Base.isempty(prop_cache::PropagationCache) = length(prop_cache) == 0

function swapterms!(prop_cache::PropagationCache)
    prop_cache.terms, prop_cache.aux_terms = prop_cache.aux_terms, prop_cache.terms
    return prop_cache
end

function zeroinactive!(prop_cache::PropagationCache)
    TT = termtype(prop_cache)
    CT = coefftype(prop_cache)

    AK.foreachindex(prop_cache.terms) do ii
        if ii > prop_cache.active_size
            prop_cache.terms[ii] = WeightedString(zero(TT), zero(CT))
        end
    end
    return prop_cache
end

function topaulisum(nq, prop_cache::PropagationCache)
    TT = termtype(prop_cache)
    CT = coefftype(prop_cache)
    # return PauliSum(nq, Dict(s.term => s.coeff for s in strings))
    psum = PauliSum(CT, nq)
    sizehint!(psum.terms, length(prop_cache.active_size))
    for wstring in viewterms(prop_cache)
        add!(psum, wstring.term, wstring.coeff)
    end
    return psum
end


function vectorpropagate(circuit, strings::Vector; kwargs...)
    return vectorpropagate(circuit, PropagationCache(strings); kwargs...).terms
end

function vectorpropagate(circuit, prop_cache::PropagationCache; min_abs_coeff=1e-10, max_weight=Inf, kwargs...)
    # assume circuit contains the parameters via freezing, or is parameter-free
    @assert countparameters(circuit) == 0 "'circuit' must be parameter-free. Consider using 'freeze()'."

    n_gates = length(circuit)
    for (i, gate) in enumerate(reverse(circuit))

        prop_cache = applygate!(gate, prop_cache)
        # println("Done with applygate!()")

        prop_cache = mergeterms!(prop_cache)
        # println("Done with mergeterms!()")
        prop_cache = truncate!(prop_cache; min_abs_coeff, max_weight)
        # println("Done with truncate!()")

        ## Thoughts:
        # - Can we merge and truncate in one?
        # - Is coeff truncation strictly worse than max cache size?
        # - Can we write a sorting function between the two sorted arrays?


    end
    # TODO: zero-out inactive parts so people are not confused
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
    AK.foreachindex(viewterms(prop_cache)) do ii
        wstring = prop_cache.terms[ii]

        term, coeff = apply(gate, wstring.term, wstring.coeff)

        prop_cache.terms[ii] = WeightedString(term, coeff)
    end

    return prop_cache
end


function PauliPropagation.getnewpaulistring(gate_mask::TT, pstr::TT) where TT
    new_pstr = PauliPropagation._bitpaulimultiply(gate_mask, pstr)

    # this counts the exponent of the imaginary unit in the new Pauli string
    im_count = PauliPropagation._calculatesignexponent(gate_mask, pstr)

    # now, instead of computing im^im_count followed by another im factor from the gate rules,
    # we do this in one step via a cheeky trick:
    sign = (im_count & 2) - 1
    # this is equivalent to sign = real( im * im^im_count)

    return new_pstr, sign
end

function applygate!(frozen_gate::FrozenGate{PauliRotation,PT}, prop_cache::PropagationCache) where {PT}

    TT = termtype(prop_cache)
    CT = coefftype(prop_cache)

    n_old = length(prop_cache)

    # unpack the frozen PauliRotation
    gate, theta = frozen_gate.gate, frozen_gate.parameter

    # this allows for faster operations
    masked_gate = PauliPropagation._tomaskedpaulirotation(gate, TT)

    # pre-compute the sine and cosine values because the are used for every Pauli string that does not commute with the gate
    cos_val = cos(theta)
    sin_val = sin(theta)

    # get the mask out because because the gate cannot be in the function when using GPU
    gate_mask = masked_gate.generator_mask
    anticommfunc = wstring -> anticommutes(wstring, gate_mask)
    flagterms!(prop_cache, anticommfunc)

    flagstoindices!(prop_cache)

    n_noncommutes = prop_cache.indices[prop_cache.active_size]

    # slit off into the same array
    n_new = n_old + n_noncommutes
    if length(prop_cache.terms) < n_new
        resize!(prop_cache, n_new)
    end

    term_view = viewterms(prop_cache)
    terms = prop_cache.terms
    flags = prop_cache.flags
    indices = prop_cache.indices
    AK.foreachindex(term_view) do ii
        # here it anticommutes
        if flags[ii]
            wstring = terms[ii]
            term = wstring.term
            coeff = wstring.coeff

            coeff1 = coeff * cos_val
            new_term, sign = PauliPropagation.getnewpaulistring(gate_mask, term)
            coeff2 = coeff * sin_val * sign

            terms[ii] = WeightedString(term, coeff1)
            terms[n_old+indices[ii]] = WeightedString(new_term, coeff2)
        end
    end
    prop_cache.active_size = n_new

    return prop_cache
end

anticommutes(wstring::WeightedString{TT,CT}, gate_generator::TT) where {TT,CT} = anticommutes(wstring.term, gate_generator)
anticommutes(term::TT, gate_generator::TT) where {TT} = !commutes(term, gate_generator)

function flagterms!(prop_cache::PropagationCache, by::F) where F
    term_view = viewterms(prop_cache)
    flags = prop_cache.flags
    AK.foreachindex(term_view) do ii
        flags[ii] = by(term_view[ii])
    end
    return
end

term(wstring::WeightedString) = wstring.term

flagstoindices!(prop_cache::PropagationCache) = flagstoindices!(prop_cache.indices, viewflags(prop_cache))
flagstoindices!(indices_dst, flags) = AK.accumulate!(+, indices_dst, flags; init=0)

function mergeterms!(prop_cache::PropagationCache)
    if isempty(prop_cache)
        return prop_cache
    end

    # Sort the vector to group identical keys together.
    AK.sort!(viewterms(prop_cache); by=term, temp=viewauxterms(prop_cache))

    # Find the start of each group in parallel.

    term_view = viewterms(prop_cache)
    flags = prop_cache.flags
    AK.foreachindex(term_view) do ii
        if ii == 1
            flags[ii] = true
        else
            flags[ii] = term_view[ii].term != term_view[ii-1].term
        end
    end

    # get the indices for where to add the coefficients
    flagstoindices!(prop_cache)
    n_unique_terms = prop_cache.indices[prop_cache.active_size]

    # early stop if all are unique 
    if n_unique_terms == prop_cache.active_size
        return prop_cache
    end

    # the reason we don't do atomic add instead is because the entry might not start off as 0
    deduplicate!(prop_cache)

    # swap terms and aux_terms
    swapterms!(prop_cache)
    prop_cache.active_size = n_unique_terms

    return prop_cache
end

function deduplicate!(prop_cache::PropagationCache)

    term_view = viewterms(prop_cache)
    aux_terms = prop_cache.aux_terms
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
            CT = typeof(term_view[ii].coeff)
            merged_coeff = zero(CT)
            for jj in ii:end_idx
                merged_coeff += term_view[jj].coeff
            end

            aux_terms[indices[ii]] = WeightedString(term_view[ii].term, merged_coeff)
        end
    end
    return
end



function truncate!(prop_cache::PropagationCache{TT,CT}; min_abs_coeff::Real, max_weight::Real=Inf) where {TT,CT}

    if isempty(prop_cache)
        return prop_cache
    end

    # TODO: can we simplify this entire function via a parallel filter!() function?

    # flag the qubits that we keep
    flagterms!(prop_cache, wstring -> (abs(tonumber(wstring.coeff)) >= min_abs_coeff) && countweight(wstring.term) <= max_weight)

    # get the new indices after deletion
    flagstoindices!(prop_cache)
    n_kept = prop_cache.indices[prop_cache.active_size]

    # pre-allocate the output array
    term_view = viewterms(prop_cache)
    aux_terms = prop_cache.aux_terms
    flags = prop_cache.flags
    indices = prop_cache.indices
    AK.foreachindex(term_view) do ii
        if flags[ii]
            aux_terms[indices[ii]] = term_view[ii]
        end
    end

    swapterms!(prop_cache)
    prop_cache.active_size = n_kept

    return prop_cache
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


function PauliPropagation.overlapwithcomputational(strings::AbstractArray{WeightedString{TT,CT}}, onebitinds) where {TT,CT}
    val = zero(CT)
    for wstring in strings
        val += tonumber(wstring.coeff) * PauliPropagation._calcsignwithones(wstring.term, onebitinds)
    end
    return val
end