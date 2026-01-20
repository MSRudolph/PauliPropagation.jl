### generics.jl
##
# This file contains the foundational functions for the `propagation` function. 
# They can be overloaded to custom gate types or custom behaviour in `specializations.jl`.
##
###
"""
    propagate(circ, pstr::PauliString, thetas=nothing; max_weight=Inf, min_abs_coeff=1e-10, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, kwargs...)

Propagate a `PauliString` through the circuit `circ` in the Heisenberg picture. 
This means that the circuit is applied to the Pauli string in reverse order, and the action of each gate is its conjugate action.
Parameters for the parametrized gates in `circ` are given by `thetas`, and need to be passed as if the circuit was applied as written in the Schrödinger picture.
If thetas are not passed, the circuit must contain only non-parametrized `StaticGates`.
Default truncations are `max_weight`, `min_abs_coeff`, `max_freq`, and `max_sins`.
`max_freq`, and `max_sins` will lead to automatic conversion if the coefficients are not already wrapped in suitable `PathProperties` objects.
A custom truncation function can be passed as `customtruncfunc` with the signature customtruncfunc(pstr::PauliStringType, coefficient)::Bool.
Further `kwargs` are passed to the lower-level functions `applymergetruncate!`, `applytoall!`, `applyandadd!`, and `apply`.
"""
function PropagationBase.propagate(circ, pstr::PauliString, thetas=nothing; max_weight=Inf, min_abs_coeff=eps(), max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, kwargs...)
    psum = PauliSum(pstr)
    return propagate(circ, psum, thetas; max_weight, min_abs_coeff, max_freq, max_sins, customtruncfunc, kwargs...)
end


function PropagationBase.propagate!(circ, pstr::PauliString, thetas=nothing; max_weight=Inf, min_abs_coeff=eps(), max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, kwargs...)
    psum = PauliSum(pstr)
    # check that max_freq and max_sins are only used a PathProperties type tracking them
    _checkfreqandsinfields(psum, max_freq, max_sins)
    return propagate(circ, psum, thetas; max_weight, min_abs_coeff, max_freq, max_sins, customtruncfunc, kwargs...)
end

"""
    propagate(circ, psum::PauliSum, thetas=nothing; max_weight=Inf, min_abs_coeff=1e-10, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, kwargs...)

Propagate a `PauliSum` through the circuit `circ` in the Heisenberg picture. 
This means that the circuit is applied to the Pauli sum in reverse order, and the action of each gate is its conjugate action.
Parameters for the parametrized gates in `circ` are given by `thetas`, and need to be passed as if the circuit was applied as written in the Schrödinger picture.
If thetas are not passed, the circuit must contain only non-parametrized `StaticGates`.
Default truncations are `max_weight`, `min_abs_coeff`, `max_freq`, and `max_sins`.
`max_freq`, and `max_sins` will lead to automatic conversion if the coefficients are not already wrapped in suitable `PathProperties` objects.
A custom truncation function can be passed as `customtruncfunc` with the signature customtruncfunc(pstr::PauliStringType, coefficient)::Bool.
Further `kwargs` are passed to the lower-level functions `applymergetruncate!`, `applytoall!`, `applyandadd!`, and `apply`.
"""
function PropagationBase.propagate(circ, psum::AbstractPauliSum, thetas=nothing; max_weight=Inf, min_abs_coeff=eps(), max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, kwargs...)
    CT = coefftype(psum)

    # if max_freq and max_sins are used, and no PathProperties used, automatically wrap the coefficients in `PauliFreqTracker` 
    psum = _check_wrapping_into_paulifreqtracker(psum, max_freq, max_sins)

    # check that max_freq and max_sins are only used a PathProperties type tracking them
    _checkfreqandsinfields(psum, max_freq, max_sins)

    # run the in-place propagation function on a deepcopy of the input psum
    psum = propagate!(circ, deepcopy(psum), thetas; max_weight, min_abs_coeff, max_freq, max_sins, customtruncfunc, kwargs...)

    # if the input psum was not a `PauliFreqTracker`, and the corresponding truncations were set,we need to unwrap the coefficients
    psum = _check_unwrap_from_paulifreqtracker(CT, psum)

    return psum
end


### TRUNCATE


function PropagationBase.truncate!(prop_cache::AbstractPauliPropagationCache; max_weight::Real=Inf, min_abs_coeff=1e-10, max_freq::Real=Inf, max_sins::Real=Inf, customtruncfunc=nothing, kwargs...)

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

