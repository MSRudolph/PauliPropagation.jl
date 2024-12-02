## This file contains the foundational functions for the `propagation` function. 
## They can be overloaded to custom gate types or custom behaviour in `specializations.jl`.

"""
    propagate(circ, pstr::PauliString, thetas; kwargs...)

Propagate a `PauliString` through the circuit `circ` in the Heisenberg picture. 
Parameters for the parametrized gates in `circ` are given by `thetas`.
`kwargs` are passed to the truncation function. Supported by default are `max_weight`, `min_abs_coeff`, `max_freq`, and `max_sins`.
A custom truncation function can be passed as `customtruncatefn` with the signature customtruncatefn(pstr::PauliStringType, coefficient)::Bool.
"""
function propagate(circ, pstr::PauliString, thetas; kwargs...)
    psum = PauliSum(pstr.nqubits, pstr)
    return propagate(circ, psum, thetas; kwargs...)
end

"""
    propagate(circ, psum::PauliSum, thetas; kwargs...)

Propagate a `PauliSum` through the circuit `circ` in the Heisenberg picture. 
Parameters for the parametrized gates in `circ` are given by `thetas`.
`kwargs` are passed to the truncation function. Supported by default are `max_weight`, `min_abs_coeff`, `max_freq`, and `max_sins`.
A custom truncation function can be passed as `customtruncatefn` with the signature customtruncatefn(pstr::PauliStringType, coefficient)::Bool.
"""
function propagate(circ, psum::PauliSum, thetas; kwargs...)
    pauli_dict = propagate!(circ, deepcopy(psum.op_dict), thetas; kwargs...)
    return PauliSum(psum.nqubits, pauli_dict)
end


"""
    propagate!(circ, psum::PauliSum, thetas; kwargs...)

Propagate a `PauliSum` through the circuit `circ` in the Heisenberg picture. 
The input `psum` will be modified.
Parameters for the parametrized gates in `circ` are given by `thetas`.
`kwargs` are passed to the truncation function. Supported by default are `max_weight`, `min_abs_coeff`, `max_freq`, and `max_sins`.
A custom truncation function can be passed as `customtruncatefn` with the signature customtruncatefn(pstr::PauliStringType, coefficient)::Bool.
"""
function propagate!(circ, psum::PauliSum, thetas; kwargs...)
    propagate!(circ, psum.op_dict, thetas; kwargs...)
    return psum
end

"""
    propagate!(circ, psum::Dict, thetas; kwargs...)

Propagate a Pauli sum of type `Dict{PauliStringType,CoeffType}` through the circuit `circ` in the Heisenberg picture. 
The input `psum` will be modified.
Parameters for the parametrized gates in `circ` are given by `thetas`.
`kwargs` are passed to the truncation function. Supported by default are `max_weight`, `min_abs_coeff`, `max_freq`, and `max_sins`.
A custom truncation function can be passed as `customtruncatefn` with the signature customtruncatefn(pstr::PauliStringType, coefficient)::Bool.
"""
function propagate!(circ, psum::Dict, thetas; kwargs...)
    # TODO: Should we have a version of this that doesn't require thetas and uses surrogate code?
    param_idx = length(thetas)

    second_psum = typeof(psum)()  # pre-allocating somehow doesn't do anything

    ## TODO:
    # - decide where to reverse the circuit
    # - verbose option  
    # - more elegant param_idx incrementation
    for gate in reverse(circ)
        psum, second_psum, param_idx = mergingapply!(gate, psum, second_psum, thetas, param_idx; kwargs...)
    end
    return psum
end

"""
    mergingapply!(gate, psum, second_psum, thetas, param_idx, args...; kwargs...)

1st-level function below `propagate!` that applies one gate to all operators in `psum`, potentially using `second_psum` in the process,
and merges everything back into `psum`. Truncations are checked here after merging.
This function can be overwritten for a custom gate if the lower-level functions `applygatetoall!`, `applygatetoone!`, and `apply` are not sufficient.
A custom truncation function can be passed as `customtruncatefn` with the signature customtruncatefn(pstr::PauliStringType, coefficient)::Bool.
"""
function mergingapply!(gate, psum, second_psum, thetas, param_idx, args...; kwargs...)

    theta = thetas[param_idx]
    psum, second_psum = applygatetoall!(gate, theta, psum, second_psum; kwargs...)

    mergeandclear!(psum, second_psum)

    checktruncationonall!(psum; kwargs...)

    if isa(gate, ParametrizedGate) && param_idx > 1  # decrement parameter index by one if it is not the last parameter
        param_idx -= 1
    end

    return psum, second_psum, param_idx
end

"""
    applygatetoall!(gate, theta psum, second_psum, args...; kwargs...)

2nd-level function below `mergingapply!` that applies one gate to all operators in `psum`, potentially using `second_psum` in the process.
This function can be overwritten for a custom gate if the lower-level functions `applygatetoone!` and `apply` are not sufficient.
"""
function applygatetoall!(gate, theta, psum, second_psum, args...; kwargs...)

    for (operator, coeff) in psum
        applygatetoone!(gate, operator, coeff, theta, psum, second_psum; kwargs...)
    end

    empty!(psum)  # empty old dict because next generation of operators should by default stored in second_psum (unless this is overwritten by a custom function)

    return second_psum, psum  # swap dicts around
end

"""
    applygatetoone!(gate, operator, coefficient, theta, psum, second_psum, args...; kwargs...)

3nd-level function below `mergingapply!` that applies one gate to one Pauli string in `psum`, potentially using `second_psum` in the process.
This function can be overwritten for a custom gate if the lower-level function `apply` is not sufficient. 
This is likely the the case if `apply` is not type-stable because it does not return a unique number of outputs. 
E.g., a Pauli gate returns 1 or 2 (operator, coefficient) outputs.
"""
@inline function applygatetoone!(gate, operator, coefficient, theta, psum, second_psum, args...; kwargs...)

    ops_and_coeffs = apply(gate, operator, theta, coefficient; kwargs...)

    for ii in 1:2:length(ops_and_coeffs)
        op, coeff = ops_and_coeffs[ii], ops_and_coeffs[ii+1]
        second_psum[ops_and_coeffs[ii]] = get(second_psum, op, 0.0) + coeff
    end

    return
end


### MERGE
"""
    mergeandclear!(psum, second_psum)

Merge `second_psum` into `psum` using the `merge` function. `merge` can be overloaded for different coefficient types.
Then clear `second_psum` for the next iteration.
"""
function mergeandclear!(psum, second_psum)
    # TODO: custom merging function beyond mergewith!
    # TODO: Potentially check for truncations at this step.
    mergewith!(merge, psum, second_psum)
    empty!(second_psum)
    return psum, second_psum
end

"""
    merge(val1, val2)

Merging two coefficients calls `+` by default unless there exists a suitable overloaded `merge` function.
"""
function merge(coeff1, coeff2)
    return coeff1 + coeff2
end


### TRUNCATE
"""
    checktruncationonall!(psum; max_weight::Real=Inf, min_abs_coeff=1e-10, max_freq::Real=Inf, max_sins::Real=Inf, kwargs...)

Check truncation conditions on all operators in `psum` and remove them if they are truncated.
This function supports the default truncations based on `max_weight`, `min_abs_coeff`, `max_freq`, and `max_sins`.
A custom truncation function can be passed as `customtruncatefn` with the signature customtruncatefn(pstr::PauliStringType, coefficient)::Bool.
"""
function checktruncationonall!(
    psum; max_weight::Real=Inf, min_abs_coeff=1e-10, max_freq::Real=Inf,
    max_sins::Real=Inf,
    kwargs...
)
    # TODO: This does currently hinder performance, even if we don't truncated
    # approx 55ms -> 66ms for the test case
    for (operator, coeff) in psum
        checktruncationonone!(
            psum, operator, coeff;
            max_weight=max_weight, min_abs_coeff=min_abs_coeff,
            max_freq=max_freq, max_sins=max_sins,
            kwargs...
        )
    end
    return
end

"""
    checktruncationonone!(
    psum, operator, coeff;
    max_weight::Real=Inf, min_abs_coeff=1e-10,
    max_freq::Real=Inf, max_sins::Real=Inf,
    customtruncatefn=nothing,
    kwargs...

Check truncation conditions one operator in `psum` and it them if it is truncated.
This function supports the default truncations based on `max_weight`, `min_abs_coeff`, `max_freq`, and `max_sins`.
A custom truncation function can be passed as `customtruncatefn` with the signature customtruncatefn(pstr::PauliStringType, coefficient)::Bool.
"""
@inline function checktruncationonone!(
    psum, operator, coeff;
    max_weight::Real=Inf, min_abs_coeff=1e-10,
    max_freq::Real=Inf, max_sins::Real=Inf,
    customtruncatefn=nothing,
    kwargs...
)
    is_truncated = false
    if truncateweight(operator, max_weight)
        is_truncated = true
    elseif truncatemincoeff(coeff, min_abs_coeff)
        is_truncated = true
    elseif truncatefrequency(coeff, max_freq)
        is_truncated = true
    elseif truncatesins(coeff, max_sins)
        is_truncated = true
    elseif !isnothing(customtruncatefn) && customtruncatefn(operator, coeff)
        is_truncated = true
    end
    if is_truncated
        delete!(psum, operator)
    end
    return
end