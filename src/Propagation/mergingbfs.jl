
function mergingbfs(circ, pstr::PauliString, thetas; kwargs...)
    psum = PauliSum(pstr.nqubits, pstr)
    return mergingbfs(circ, psum, thetas; kwargs...)
end

function mergingbfs!(circ, pstr::PauliString, thetas; kwargs...)
    return mergingbfs(circ, pstr, thetas; kwargs...)
end

function mergingbfs(circ, psum::PauliSum, thetas; kwargs...)
    pauli_dict = mergingbfs!(circ, deepcopy(psum.op_dict), thetas; kwargs...)
    return PauliSum(psum.nqubits, pauli_dict)
end

function mergingbfs!(circ, psum::PauliSum, thetas; kwargs...)
    pauli_dict = mergingbfs!(circ, psum.op_dict, thetas; kwargs...)
    return PauliSum(psum.nqubits, pauli_dict)
end


function mergingbfs!(circ, d, thetas; kwargs...)
    # TODO: Should we have a version of this that doesn't require thetas and uses surrogate code?
    param_idx = length(thetas)

    second_d = typeof(d)()  # pre-allocating somehow doesn't do anything

    ## TODO:
    # - decide where to reverse the circuit
    # - verbose option  
    for gate in reverse(circ)
        d, second_d, param_idx = mergingapply(gate, d, second_d, thetas, param_idx; kwargs...)
    end
    return d
end


function mergingapply(gate, operator_dict::Dict, new_operator_dict::Dict, thetas, param_idx, args...; kwargs...)

    # param_idx is decremented by one if the gate is a Pauli gate
    operator_dict, new_operator_dict, param_idx = applygatetoall!(gate, thetas, param_idx, operator_dict, new_operator_dict)

    mergeandclear!(operator_dict, new_operator_dict)

    checktruncationonall!(operator_dict; kwargs...)

    return operator_dict, new_operator_dict, param_idx
end

### APPLY

function applygatetoall!(gate, thetas, param_idx, operator_dict, new_operator_dict, args...; kwargs...)
    theta = thetas[param_idx]
    for (operator, coeff) in operator_dict
        applygatetoone!(gate, operator, coeff, theta, param_idx, operator_dict, new_operator_dict; kwargs...)
    end

    if isa(gate, ParametrizedGate)  # decrement parameter index by one
        param_idx -= 1
    end

    empty!(operator_dict)  # empty old dict because next generation of operators should by default stored in new_operator_dict (unless this is a overwritten by a custom function)

    return new_operator_dict, operator_dict, param_idx  # swap dicts around
end

@inline function applygatetoone!(gate, operator, coefficient, theta, param_idx, operator_dict, new_operator_dict, args...; kwargs...)

    ops_and_coeffs = apply(gate, operator, theta, coefficient)

    for ii in 1:2:length(ops_and_coeffs)
        op, coeff = ops_and_coeffs[ii], ops_and_coeffs[ii+1]
        new_operator_dict[ops_and_coeffs[ii]] = get(new_operator_dict, op, 0.0) + coeff
    end

    return
end

### PAULI GATES

function applygatetoall!(gate::PauliGateUnion, thetas, param_idx, operator_dict, new_operator_dict, args...; kwargs...)  # TODO: there is a lot of code duplication. Can we write a more general function?
    theta = thetas[param_idx]
    for (operator, coeff) in operator_dict
        applygatetoone!(gate, operator, coeff, theta, param_idx, operator_dict, new_operator_dict; kwargs...)
    end

    param_idx -= 1

    return operator_dict, new_operator_dict, param_idx
end

@inline function applygatetoone!(gate::PauliGateUnion, operator, coefficient, theta, param_idx, operator_dict, new_operator_dict, args...; kwargs...)

    if commutes(gate, operator)
        return
    end

    operator, coeff1, new_oper, coeff2 = applynoncummuting(gate, operator, theta, coefficient; param_idx=param_idx)  # TODO: remove the param_idx from here because it is only specific to the Surrogate

    operator_dict[operator] = coeff1
    new_operator_dict[new_oper] = coeff2

    return
end


### Amplitude Damping Noise

function applygatetoall!(gate::AmplitudeDampingNoise, thetas, param_idx, operator_dict, new_operator_dict, args...; kwargs...)  # TODO: there is a lot of code duplication. Can we write a more general function?
    theta = thetas[param_idx]
    for (operator, coeff) in operator_dict
        applygatetoone!(gate, operator, coeff, theta, param_idx, operator_dict, new_operator_dict; kwargs...)
    end

    param_idx -= 1

    return operator_dict, new_operator_dict, param_idx
end

@inline function applygatetoone!(gate::AmplitudeDampingNoise, operator, coefficient, theta, param_idx, operator_dict, new_operator_dict, args...; kwargs...)

    if actsdiagonally(gate, operator)
        operator, coeff = diagonalapply(gate, operator, theta, coefficient)
        operator_dict[operator] = coeff
        return
    end

    operator, coeff1, new_oper, coeff2 = splitapply(gate, operator, theta, coefficient)

    operator_dict[operator] = coeff1
    new_operator_dict[new_oper] = coeff2

    return
end

### Clifford gates

## NOTE: Currently, a slightly more optimized applygatetoall! function that gets the map_array once and passes that to applygatetoone! is not significantly faster than the general function


### MERGE

# TODO: custom merging function beyond mergewith!
function mergeandclear!(operator_dict::Dict, new_operator_dict::Dict)
    mergewith!(merge, operator_dict, new_operator_dict)
    empty!(new_operator_dict)
    return operator_dict, new_operator_dict
end

function merge(val1, val2)
    return val1 + val2
end

function merge(pth1::PathProperties, pth2::PathProperties)
    pth1.coeff = merge(pth1.coeff, pth2.coeff)
    pth1.ncos = min(pth1.ncos, pth2.ncos)
    pth1.nsins = min(pth1.nsins, pth2.nsins)
    pth1.freq = min(pth1.freq, pth2.freq)
    return pth1
end

function merge(coeff1::Number, coeff2::Number)
    return coeff1 + coeff2
end

### TRUNCATE

function checktruncationonall!(
    operator_dict; max_weight::Real=Inf, min_abs_coeff=0.0, max_freq::Real=Inf,
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

@inline function checktruncationonone!(
    operator_dict, operator, coeff;
    max_weight::Real=Inf, min_abs_coeff=0.0,
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