"""
    mergingbfs(circ, pstr::PauliString, thetas; kwargs...)

Perform merging breadth-first search (BFS) simulation of a `PauliString` propagating through the circuit `circ` in the Heisenberg picture. 
Parameters for the parametrized gates in `circ` are given by `thetas`.
`kwargs` are passed to the truncation function. Supported by default are `max_weight`, `min_abs_coeff`, `max_freq`, and `max_sins`.
A custom truncation function can be passed as `customtruncatefn` with the signature customtruncatefn(pstr::PauliStringType, coefficient)::Bool.
"""
function mergingbfs(circ, pstr::PauliString, thetas; kwargs...)
    psum = PauliSum(pstr.nqubits, pstr)
    return mergingbfs(circ, psum, thetas; kwargs...)
end

"""
    mergingbfs(circ, psum::PauliSum, thetas; kwargs...)

Perform merging breadth-first search (BFS) simulation of a `PauliSum` propagating through the circuit `circ` in the Heisenberg picture. 
Parameters for the parametrized gates in `circ` are given by `thetas`.
`kwargs` are passed to the truncation function. Supported by default are `max_weight`, `min_abs_coeff`, `max_freq`, and `max_sins`.
A custom truncation function can be passed as `customtruncatefn` with the signature customtruncatefn(pstr::PauliStringType, coefficient)::Bool.
"""
function mergingbfs(circ, psum::PauliSum, thetas; kwargs...)
    pauli_dict = mergingbfs!(circ, deepcopy(psum.op_dict), thetas; kwargs...)
    return PauliSum(psum.nqubits, pauli_dict)
end


"""
    mergingbfs!(circ, psum::PauliSum, thetas; kwargs...)

Perform in-place merging breadth-first search (BFS) simulation of a `PauliSum` propagating through the circuit `circ` in the Heisenberg picture. 
The input `psum` will be modified.
Parameters for the parametrized gates in `circ` are given by `thetas`.
`kwargs` are passed to the truncation function. Supported by default are `max_weight`, `min_abs_coeff`, `max_freq`, and `max_sins`.
A custom truncation function can be passed as `customtruncatefn` with the signature customtruncatefn(pstr::PauliStringType, coefficient)::Bool.
"""
function mergingbfs!(circ, psum::PauliSum, thetas; kwargs...)
    pauli_dict = mergingbfs!(circ, psum.op_dict, thetas; kwargs...)
    return PauliSum(psum.nqubits, pauli_dict)
end

"""
    mergingbfs!(circ, d::Dict, thetas; kwargs...)

Perform in-place merging breadth-first search (BFS) simulation of a `Dict{PauliStringType,CoeffType}` propagating through the circuit `circ` in the Heisenberg picture. 
The input `psum` will be modified.
Parameters for the parametrized gates in `circ` are given by `thetas`.
`kwargs` are passed to the truncation function. Supported by default are `max_weight`, `min_abs_coeff`, `max_freq`, and `max_sins`.
A custom truncation function can be passed as `customtruncatefn` with the signature customtruncatefn(pstr::PauliStringType, coefficient)::Bool.
"""
function mergingbfs!(circ, d::Dict, thetas; kwargs...)
    # TODO: Should we have a version of this that doesn't require thetas and uses surrogate code?
    param_idx = length(thetas)

    second_d = typeof(d)()  # pre-allocating somehow doesn't do anything

    ## TODO:
    # - decide where to reverse the circuit
    # - verbose option  
    # - more elegant param_idx incrementation
    for gate in reverse(circ)
        d, second_d, param_idx = mergingapply(gate, d, second_d, thetas, param_idx; kwargs...)
    end
    return d
end

"""
    mergingapply(gate, operator_dict::Dict, new_operator_dict::Dict, thetas, param_idx, args...; kwargs...)

1st-level function below `mergingbfs!` that applies one gate to all operators in `operator_dict`, potentially using `new_operator_dict` in the process,
and merges everything back into `operator_dict`. Truncations are checked here after merging.
This function can be overwritten for a custom gate if the lower-level functions `applygatetoall!`, `applygatetoone!`, and `apply` are not sufficient.
A custom truncation function can be passed as `customtruncatefn` with the signature customtruncatefn(pstr::PauliStringType, coefficient)::Bool.
"""
function mergingapply(gate, operator_dict, new_operator_dict, thetas, param_idx, args...; kwargs...)

    theta = thetas[param_idx]
    operator_dict, new_operator_dict = applygatetoall!(gate, theta, operator_dict, new_operator_dict; kwargs...)

    mergeandclear!(operator_dict, new_operator_dict)

    checktruncationonall!(operator_dict; kwargs...)

    if isa(gate, ParametrizedGate) && param_idx > 1  # decrement parameter index by one if it is not the last parameter
        param_idx -= 1
    end

    return operator_dict, new_operator_dict, param_idx
end

"""
    applygatetoall!(gate, theta operator_dict, new_operator_dict, args...; kwargs...)

2nd-level function below `mergingapply!` that applies one gate to all operators in `operator_dict`, potentially using `new_operator_dict` in the process.
This function can be overwritten for a custom gate if the lower-level functions `applygatetoone!` and `apply` are not sufficient.
"""
function applygatetoall!(gate, theta, operator_dict, new_operator_dict, args...; kwargs...)

    for (operator, coeff) in operator_dict
        applygatetoone!(gate, operator, coeff, theta, operator_dict, new_operator_dict; kwargs...)
    end

    empty!(operator_dict)  # empty old dict because next generation of operators should by default stored in new_operator_dict (unless this is a overwritten by a custom function)

    return new_operator_dict, operator_dict  # swap dicts around
end

"""
    applygatetoone!(gate, operator, coefficient, theta, operator_dict, new_operator_dict, args...; kwargs...)

3nd-level function below `mergingapply!` that applies one gate to one Pauli string in `operator_dict`, potentially using `new_operator_dict` in the process.
This function can be overwritten for a custom gate if the lower-level function `apply` is not sufficient. 
This is likely the the case if `apply` is not type-stable because it does not return a unique number of outputs. E.g., a Pauli gate returns 1 or 2 (operator, coefficient) outputs.
"""
@inline function applygatetoone!(gate, operator, coefficient, theta, operator_dict, new_operator_dict, args...; kwargs...)

    ops_and_coeffs = apply(gate, operator, theta, coefficient; kwargs...)

    for ii in 1:2:length(ops_and_coeffs)
        op, coeff = ops_and_coeffs[ii], ops_and_coeffs[ii+1]
        new_operator_dict[ops_and_coeffs[ii]] = get(new_operator_dict, op, 0.0) + coeff
    end

    return
end

### PAULI GATES
"""
    applygatetoall!(gate::PauliGateUnion, thetas, operator_dict, new_operator_dict, args...; kwargs...)

Overload of `applygatetoall!` for `PauliGate` and `FastPauliGate` gates. Both `operator_dict` and `new_operator_dict` contain operators which will be merged later.
"""
function applygatetoall!(gate::PauliGateUnion, theta, operator_dict, new_operator_dict, args...; kwargs...)
    # TODO: there is a lot of code duplication. Can we write a more general function?

    for (operator, coeff) in operator_dict
        applygatetoone!(gate, operator, coeff, theta, operator_dict, new_operator_dict; kwargs...)
    end

    return operator_dict, new_operator_dict
end

"""
    applygatetoone!(gate::PauliGateUnion, operator, coefficient, theta, operator_dict, new_operator_dict, args...; kwargs...)

Overload of `applygatetoone!` for `PauliGate` and `FastPauliGate` gates. Checks for commutation of `gate` and `operator`, and applies the gate to the operator if they don't.
"""
@inline function applygatetoone!(gate::PauliGateUnion, operator, coefficient, theta, operator_dict, new_operator_dict, args...; kwargs...)

    if commutes(gate, operator)
        return
    end

    operator, coeff1, new_oper, coeff2 = applynoncummuting(gate, operator, theta, coefficient; kwargs...)

    operator_dict[operator] = coeff1
    new_operator_dict[new_oper] = coeff2

    return
end


### Amplitude Damping Noise
"""
    applygatetoall!(gate::PauliGateUnion, thetas, operator_dict, new_operator_dict, args...; kwargs...)

Overload of `applygatetoall!` for `AmplitudeDampingNoise` gates. Both `operator_dict` and `new_operator_dict` contain operators which will be merged later.
"""
function applygatetoall!(gate::AmplitudeDampingNoise, theta, operator_dict, new_operator_dict, args...; kwargs...)
    # TODO: there is a lot of code duplication. Can we write a more general function? 

    for (operator, coeff) in operator_dict
        applygatetoone!(gate, operator, coeff, theta, operator_dict, new_operator_dict; kwargs...)
    end

    return operator_dict, new_operator_dict
end

"""
    applygatetoone!(gate::PauliGateUnion, operator, coefficient, theta, operator_dict, new_operator_dict, args...; kwargs...)

Overload of `applygatetoone!` for `AmplitudeDampingNoise` gates. Checks for whether `gate` will cause splitting and has tailored logic.
"""
@inline function applygatetoone!(gate::AmplitudeDampingNoise, operator, coefficient, theta, operator_dict, new_operator_dict, args...; kwargs...)

    if actsdiagonally(gate, operator)
        operator, coeff = diagonalapply(gate, operator, theta, coefficient; kwargs...)
        operator_dict[operator] = coeff
        return
    end

    operator, coeff1, new_oper, coeff2 = splitapply(gate, operator, theta, coefficient; kwargs...)

    operator_dict[operator] = coeff1
    new_operator_dict[new_oper] = coeff2

    return
end

### Clifford gates

"""
    applygatetoone!(gate::CliffordGate, operator, coefficient, theta, operator_dict, new_operator_dict, args...; kwargs...)

Overload of `applygatetoone!` for `CliffordGate` gates. Simplified logic for readability.
"""
@inline function applygatetoone!(gate::CliffordGate, operator, coefficient, theta, operator_dict, new_operator_dict, args...; kwargs...)

    op, coeff = apply(gate, operator, theta, coefficient; kwargs...)
    new_operator_dict[op] = coeff

    return
end

### MERGE
"""
    mergeandclear!(operator_dict, new_operator_dict)

Merge `new_operator_dict` into `operator_dict` using the `merge` function. `merge` can be overloaded for different coefficient types.
Then clear `new_operator_dict` for the next iteration.
"""
function mergeandclear!(operator_dict::Dict, new_operator_dict::Dict)
    # TODO: custom merging function beyond mergewith!
    # TODO: Potentially check for truncations at this step.
    mergewith!(merge, operator_dict, new_operator_dict)
    empty!(new_operator_dict)
    return operator_dict, new_operator_dict
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
    checktruncationonall!(operator_dict; max_weight::Real=Inf, min_abs_coeff=1e-10, max_freq::Real=Inf, max_sins::Real=Inf, kwargs...)

Check truncation conditions on all operators in `operator_dict` and remove them if they are truncated.
This function supports the default truncations based on `max_weight`, `min_abs_coeff`, `max_freq`, and `max_sins`.
A custom truncation function can be passed as `customtruncatefn` with the signature customtruncatefn(pstr::PauliStringType, coefficient)::Bool.
"""
function checktruncationonall!(
    operator_dict; max_weight::Real=Inf, min_abs_coeff=1e-10, max_freq::Real=Inf,
    max_sins::Real=Inf,
    kwargs...
)
    # TODO: This does currently hinder performance, even if we don't truncated
    # approx 55ms -> 66ms for the test case
    for (operator, coeff) in operator_dict
        checktruncationonone!(
            operator_dict, operator, coeff;
            max_weight=max_weight, min_abs_coeff=min_abs_coeff,
            max_freq=max_freq, max_sins=max_sins,
            kwargs...
        )
    end
    return
end

"""
    checktruncationonone!(
    operator_dict, operator, coeff;
    max_weight::Real=Inf, min_abs_coeff=1e-10,
    max_freq::Real=Inf, max_sins::Real=Inf,
    customtruncatefn=nothing,
    kwargs...

Check truncation conditions one operator in `operator_dict` and it them if it is truncated.
This function supports the default truncations based on `max_weight`, `min_abs_coeff`, `max_freq`, and `max_sins`.
A custom truncation function can be passed as `customtruncatefn` with the signature customtruncatefn(pstr::PauliStringType, coefficient)::Bool.
"""
@inline function checktruncationonone!(
    operator_dict, operator, coeff;
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
        delete!(operator_dict, operator)
    end
    return
end