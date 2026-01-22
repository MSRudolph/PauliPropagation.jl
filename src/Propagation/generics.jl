### generics.jl
##
# This file contains the foundational functions for the `propagation` function. 
# They can be overloaded to custom gate types or custom behaviour in `specializations.jl`.
##
###
"""
    propagate(circ, pstr::PauliString, thetas=nothing; min_abs_coeff=eps(), max_weight=Inf, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, heisenberg=true, kwargs...)

Propagate a `PauliString` through the circuit `circ` in the Heisenberg picture. 
This means that the circuit is applied to the Pauli string in reverse order, and the action of each gate is its conjugate action.
Parameters for the parametrized gates in `circ` are given by `thetas`, and need to be passed as if the circuit was applied as written in the Schrödinger picture.
If thetas are not passed, the circuit must contain only non-parametrized `StaticGates`.
Default truncations are `min_abs_coeff`, `max_weight`, `max_freq`, and `max_sins`.
`max_freq`, and `max_sins` will lead to automatic conversion if the coefficients are not already wrapped in suitable `PathProperties` objects.
A custom truncation function can be passed as `customtruncfunc` with the signature customtruncfunc(pstr::PauliStringType, coefficient)::Bool.
Further `kwargs` are passed to the lower-level functions `applymergetruncate!`, `applytoall!`, and `apply`.
"""
function propagate(circuit, pstr::PauliString, thetas=nothing; min_abs_coeff=eps(), max_weight=Inf, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, heisenberg=true, kwargs...)
    psum = PauliSum(pstr)
    return propagate(circuit, psum, thetas; min_abs_coeff, max_weight, max_freq, max_sins, customtruncfunc, kwargs...)
end

"""
    propagate!(circ, pstr::PauliString, thetas=nothing; min_abs_coeff=eps(), max_weight=Inf, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, heisenberg=true, kwargs...)

In-place version of `propagate()` for a `PauliString`.
This is only a convenience function, because the `PauliString` is converted into a `PauliSum` internally.
If `max_freq`, and `max_sins` are used without the coefficients already being wrapped in suitable `PathProperties` objects, an error is thrown.
"""
function propagate!(circuit, pstr::PauliString, thetas=nothing; max_weight=Inf, min_abs_coeff=eps(), max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, heisenberg=true, kwargs...)
    psum = PauliSum(pstr)
    # check that max_freq and max_sins are only used a PathProperties type tracking them
    _checkfreqandsinfields(psum, max_freq, max_sins)
    return propagate(circuit, psum, thetas; min_abs_coeff, max_weight, max_freq, max_sins, customtruncfunc, kwargs...)
end

"""
    propagate(circuit, psum::AbstractPauliSum, thetas=nothing; min_abs_coeff=eps(), max_weight=Inf, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, heisenberg=true, kwargs...)

Propagate a Pauli sum `psum` through the circuit `circ` in the Heisenberg picture. 
This means that the circuit is applied to the Pauli sum in reverse order, and the action of each gate is its conjugate action.
The Pauli sum `psum` is deepcopied and passed into the in-place propagation function `propagate!()`.
Parameters for the parametrized gates in `circ` are given by `thetas`, and need to be passed as if the circuit was applied as written in the Schrödinger picture.
If thetas are not passed, the circuit must contain only non-parametrized `StaticGates`.
Default truncations are `min_abs_coeff`, `max_weight`, `max_freq`, and `max_sins`.
`max_freq`, and `max_sins` will lead to automatic conversion if the coefficients are not already wrapped in suitable `PathProperties` objects.
A custom truncation function can be passed as `customtruncfunc` with the signature customtruncfunc(pstr::PauliStringType, coefficient)::Bool.
Further `kwargs` are passed to the lower-level functions `applymergetruncate!`, `applytoall!`, and `apply`.
"""
function propagate(circuit, psum::AbstractPauliSum, thetas=nothing; max_weight=Inf, min_abs_coeff=eps(), max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, heisenberg=true, kwargs...)
    CT = coefftype(psum)

    # if max_freq and max_sins are used, and no PathProperties used, automatically wrap the coefficients in `PauliFreqTracker` 
    psum = _check_wrapping_into_paulifreqtracker(psum, max_freq, max_sins)

    # check that max_freq and max_sins are only used a PathProperties type tracking them
    _checkfreqandsinfields(psum, max_freq, max_sins)

    # run the in-place propagation function on a deepcopy of the input psum
    psum = propagate!(circuit, deepcopy(psum), thetas; max_weight, min_abs_coeff, max_freq, max_sins, customtruncfunc, kwargs...)

    # if the input psum was not a `PauliFreqTracker`, and the corresponding truncations were set,we need to unwrap the coefficients
    psum = _check_unwrap_from_paulifreqtracker(CT, psum)

    return psum
end

"""
    propagate!(circuit, psum::AbstractPauliSum, thetas=nothing; min_abs_coeff=eps(), max_weight=Inf, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, heisenberg=true, kwargs...)

In-place version of `propagate()` for an `AbstractPauliSum`.
"""
function propagate!(circuit, term_sum::PS, params=nothing; max_weight=Inf, min_abs_coeff=eps(), max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, heisenberg=true, kwargs...) where PS<:AbstractPauliSum
    prop_cache = propagate!(circuit, PropagationCache(term_sum), params; max_weight=Inf, min_abs_coeff=eps(), max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, kwargs...)
    # default constructors from a cache back to our Pauli sums are defined 
    return PS(prop_cache)
end

"""
    propagate!(circuit, prop_cache::AbstractPauliPropagationCache, thetas=nothing; min_abs_coeff=eps(), max_weight=Inf, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, heisenberg=true, kwargs...)

In-place propagation of an `AbstractPauliPropagationCache` through the circuit `circ` in the Heisenberg picture.
"""
function propagate!(circuit, prop_cache::AbstractPauliPropagationCache, thetas=nothing; max_weight=Inf, min_abs_coeff=eps(), max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, heisenberg=true, kwargs...)

    # if circuit is actually a single gate, promote it to a list [gate]
    # similarly the params if it is a single number
    circuit, params = PropagationBase._promotecircandparams(circuit, params)

    # if params is nothing, the circuit must contain only StaticGates
    # also check if the length of params equals the number of parametrized gates
    PropagationBase._checkcircandparams(circuit, params)

    if heisenberg
        circuit, params = reverse(circuit), reverse(params)
    else
        circuit, params = toschrodinger(circuit, params)
    end

    return PropagationBase.propagate!(circuit, prop_cache, thetas; max_weight=Inf, min_abs_coeff=eps(), max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, kwargs...)
end

function PropagationBase.applymergetruncate!(gate, prop_cache::AbstractPauliPropagationCache, args...; kwargs...)
    # overload inside PauliPropagation.jl to skip merging for some gates
    # irrelevant for the Dict-based PauliSum, because there the main and aux sums are already swapped
    # but skipping merging is great for the VectorPauliSum which is limited by sorting

    applytoall!(gate, prop_cache, args...; kwargs...)

    if _requiresmerging(gate)
        merge!(prop_cache; kwargs...)
    end

    truncate!(prop_cache; kwargs...)

    return
end

_requiresmerging(gate) = true
_requiresmerging(::CliffordGate) = false
_requiresmerging(::PauliNoise) = false

### TRUNCATE

"""
truncate!(psum::AbstractPauliSum; min_abs_coeff=eps(), max_weight=Inf, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, kwargs...)    
truncate!(prop_cache::AbstractPauliPropagationCache; min_abs_coeff=eps(), max_weight=Inf, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, kwargs...)

Truncation function for `AbstractPauliPropagationCache`s that combines multiple truncation criteria.
The default truncation criteria are:
- `min_abs_coeff`: Truncates Pauli strings with absolute coefficient below this value.
- `max_weight`: Truncates Pauli strings with weight (number of non-identity Paulis) above this value.
- `max_freq`: Truncates Pauli strings with frequency (number of cosine factors in coefficient) above this value.
- `max_sins`: Truncates Pauli strings with number of sine factors in coefficient above this value.
A custom truncation function can be passed as `customtruncfunc` with the signature customtruncfunc(pstr, coeff)::Bool.

This function combines all truncation criteria into a single truncation function `truncfunc()` calls PropagationBase.truncate!(truncfunc, prop_cache).
"""
function PropagationBase.truncate!(prop_cache::AbstractPauliPropagationCache; min_abs_coeff::Real=eps(), max_weight::Real=Inf, max_freq::Real=Inf, max_sins::Real=Inf, customtruncfunc=nothing, kwargs...)

    function truncfunc(pstr, coeff)
        is_truncated = false
        if truncateweight(pstr, max_weight)
            is_truncated = true
        elseif truncatemincoeff(coeff, min_abs_coeff)
            is_truncated = true
        elseif truncatefrequency(coeff, max_freq)
            is_truncated = true
        elseif truncatesins(coeff, max_sins)
            is_truncated = true
        elseif !isnothing(customtruncfunc) && customtruncfunc(pstr, coeff)
            is_truncated = true
        end

        return is_truncated
    end

    prop_cache = truncate!(truncfunc, prop_cache)

    return prop_cache
end

function PropagationBase.truncate!(psum::AbstractPauliSum; min_abs_coeff::Real=eps(), max_weight::Real=Inf, max_freq::Real=Inf, max_sins::Real=Inf, customtruncfunc=nothing, kwargs...)

    function truncfunc(pstr, coeff)
        is_truncated = false
        if truncateweight(pstr, max_weight)
            is_truncated = true
        elseif truncatemincoeff(coeff, min_abs_coeff)
            is_truncated = true
        elseif truncatefrequency(coeff, max_freq)
            is_truncated = true
        elseif truncatesins(coeff, max_sins)
            is_truncated = true
        elseif !isnothing(customtruncfunc) && customtruncfunc(pstr, coeff)
            is_truncated = true
        end

        return is_truncated
    end

    term_sum = truncate!(truncfunc, psum; kwargs...)
    return term_sum
end
