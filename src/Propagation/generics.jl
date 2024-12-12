### generics.jl
##
# This file contains the foundational functions for the `propagation` function. 
# They can be overloaded to custom gate types or custom behaviour in `specializations.jl`.
##
###
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
    pauli_dict = propagate!(circ, deepcopy(psum.terms), thetas; kwargs...)
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
    propagate!(circ, psum.terms, thetas; kwargs...)
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

    aux_psum = typeof(psum)()  # pre-allocating somehow doesn't do anything

    ## TODO:
    # - decide where to reverse the circuit
    # - verbose option  
    # - more elegant param_idx incrementation
    for gate in reverse(circ)
        psum, aux_psum, param_idx = applymergetruncate!(gate, psum, aux_psum, thetas, param_idx; kwargs...)
    end
    return psum
end

# TODO: somehow propagate is not type stable

"""
    applymergetruncate!(gate, psum, aux_psum, thetas, param_idx, args...; kwargs...)

1st-level function below `propagate!` that applies one gate to all Pauli strings in `psum`, potentially using `aux_psum` in the process,
and merges everything back into `psum`. Truncations are checked here after merging.
This function can be overwritten for a custom gate if the lower-level functions `applygatetoall!`, `applygatetoone!`, and `apply` are not sufficient.
A custom truncation function can be passed as `customtruncatefn` with the signature customtruncatefn(pstr::PauliStringType, coefficient)::Bool.
"""
function applymergetruncate!(gate, psum, aux_psum, thetas, param_idx, args...; kwargs...)

    # Pick out the theta of the next parametrized gate. Will not be used by static gates, even though it is passed.
    theta = thetas[param_idx]

    # Apply the gate to all Pauli strings in psum, potentially writing into auxillary aux_psum in the process.
    # The pauli sums will be changed in-place
    applytoall!(gate, theta, psum, aux_psum; kwargs...)

    # Any contents of psum and aux_psum are merged into the larger of the two, which is returned as psum.
    # The other is cleared and returned as aux_psum.
    psum, aux_psum = mergeandclear!(psum, aux_psum)

    # Check truncation conditions on all Pauli strings in psum and remove them if they are truncated.
    psum = checktruncationonall!(psum; kwargs...)

    # If the gate was parametrized and used the theta, decrement theta index by one.
    # Don't go below index 1 because a theta will still be picked out for all remaining static gates.
    if isa(gate, ParametrizedGate) && param_idx > 1
        param_idx -= 1
    end

    return psum, aux_psum, param_idx
end

"""
    applytoall!(gate, theta psum, aux_psum, args...; kwargs...)

2nd-level function below `applymergetruncate!` that applies one gate to all Pauli strings in `psum`, moving results into `aux_psum` in the process.
This function can be overwritten for a custom gate if the lower-level functions `applygatetoone!` and `apply` are not sufficient.
In particular, this function can be used to manipulate both `psum` and `aux_psum` at the same time to reduce memory movement.
Note that manipulating `psum` on anything other than the current Pauli string will likely lead to errors.
"""
function applytoall!(gate, theta, psum, aux_psum, args...; kwargs...)

    # Loop over all Pauli strings in psum and apply the gate to them.
    for (pstr, coeff) in psum
        # apply gate to one Pauli string, move new Pauli strings to aux_psum
        applyandadd!(gate, pstr, coeff, theta, aux_psum; kwargs...)
    end
    # Empty psum because everything was moved into aux_psum. They will later be swapped.
    # If we want to reduce unnecessary Pauli string movement, we can overload applygatetoall!()
    empty!(psum)

    return
end

"""
    applyandadd!(gate, pstr, coefficient, theta, aux_psum, args...; kwargs...)

3nd-level function below `applymergetruncate!` that applies one gate to one Pauli string in `psum`, moving results into `aux_psum` in the process.
This function can be overwritten for a custom gate if the lower-level function `apply` is not sufficient. 
This is likely the the case if `apply` is not type-stable because it does not return a unique number of outputs. 
E.g., a Pauli gate returns 1 or 2 (pstr, coefficient) outputs.
"""
@inline function applyandadd!(gate, pstr, coefficient, theta, aux_psum, args...; kwargs...)

    # Get the (potentially new) pauli strings and their coefficients like (pstr1, coeff1, pstr2, coeff2, ...)
    pstrs_and_coeffs = apply(gate, pstr, theta, coefficient; kwargs...)

    for ii in 1:2:length(pstrs_and_coeffs)
        # Itererate over the pairs of pstr and coeff
        new_pstr, new_coeff = pstrs_and_coeffs[ii], pstrs_and_coeffs[ii+1]
        # Store the new_pstr and coeff in the aux_psum, add to existing coeff if new_pstr already exists there
        add!(aux_psum, new_pstr, new_coeff; precision=0)
    end

    return
end


### MERGE
"""
    mergeandclear!(psum, aux_psum)

Merge `aux_psum` into `psum` using the `merge` function. `merge` can be overloaded for different coefficient types.
Then clear `aux_psum` for the next iteration.
"""
function mergeandclear!(psum, aux_psum)
    # merge the smaller dict into the larger one
    if length(psum) < length(aux_psum)
        psum, aux_psum = aux_psum, psum
    end
    # TODO: custom merging function beyond mergewith!
    # TODO: Potentially check for truncations at this step.
    mergewith!(merge, psum, aux_psum)
    empty!(aux_psum)
    return psum, aux_psum
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

Check truncation conditions on all Pauli strings in `psum` and remove them if they are truncated.
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
    for (pstr, coeff) in psum
        checktruncationonone!(
            psum, pstr, coeff;
            max_weight=max_weight, min_abs_coeff=min_abs_coeff,
            max_freq=max_freq, max_sins=max_sins,
            kwargs...
        )
    end
    return psum
end

"""
    checktruncationonone!(
    psum, pstr, coeff;
    max_weight::Real=Inf, min_abs_coeff=1e-10,
    max_freq::Real=Inf, max_sins::Real=Inf,
    customtruncatefn=nothing,
    kwargs...

Check truncation conditions one Pauli string in `psum` and it them if it is truncated.
This function supports the default truncations based on `max_weight`, `min_abs_coeff`, `max_freq`, and `max_sins`.
A custom truncation function can be passed as `customtruncatefn` with the signature customtruncatefn(pstr::PauliStringType, coefficient)::Bool.
"""
@inline function checktruncationonone!(
    psum, pstr, coeff;
    max_weight::Real=Inf, min_abs_coeff=1e-10,
    max_freq::Real=Inf, max_sins::Real=Inf,
    customtruncatefn=nothing,
    kwargs...
)
    is_truncated = false
    if truncateweight(pstr, max_weight)
        is_truncated = true
    elseif truncatemincoeff(coeff, min_abs_coeff)
        is_truncated = true
    elseif truncatefrequency(coeff, max_freq)
        is_truncated = true
    elseif truncatesins(coeff, max_sins)
        is_truncated = true
    elseif !isnothing(customtruncatefn) && customtruncatefn(pstr, coeff)
        is_truncated = true
    end
    if is_truncated
        delete!(psum, pstr)
    end
    return
end