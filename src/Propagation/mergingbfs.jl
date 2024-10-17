
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


function mergingapply(gate, operator_dict::Dict, new_operator_dict::Dict, thetas, param_idx, args...; customtruncationfunction=nothing, min_abs_coeff=0.0, kwargs...)

    # param_idx is decremented by one if the gate is a Pauli gate
    operator_dict, new_operator_dict, param_idx = applygatetoall!(gate, thetas, param_idx, operator_dict, new_operator_dict)

    mergeandclear!(operator_dict, new_operator_dict)

    checktruncationonall!(operator_dict; kwargs...)

    if !isnothing(customtruncationfunction)
        customtruncationfunction(operator_dict, param_idx)  # changes in-place
    end

    return operator_dict, new_operator_dict, param_idx
end

### APPLY

function applygatetoall!(gate, thetas, param_idx, operator_dict, new_operator_dict, args...; kwargs...)
    theta = 0.0  # This gate seems to not use parameters
    for (operator, coeff) in operator_dict
        applygatetoone!(gate, operator, coeff, theta, param_idx, operator_dict, new_operator_dict; kwargs...)  #  max_weight=here_max_weight, 
    end
    return operator_dict, new_operator_dict, param_idx
end

function applygatetoall!(gate::PauliGateUnion, thetas, param_idx, operator_dict, new_operator_dict, args...; kwargs...)
    theta = thetas[param_idx]
    for (operator, coeff) in operator_dict
        applygatetoone!(gate, operator, coeff, theta, param_idx, operator_dict, new_operator_dict; kwargs...)  #  max_weight=here_max_weight, 
    end
    return operator_dict, new_operator_dict, param_idx - 1  # decrement parameter index by one
end

function applygatetoall!(gate::CliffordGate, thetas, param_idx, operator_dict, new_operator_dict, args...; max_weight::Real=Inf, kwargs...)

    map_array = default_clifford_map[gate.symbol]

    for (oper, coeff) in operator_dict
        applygatetoone!(gate, oper, coeff, thetas, param_idx, map_array, operator_dict, new_operator_dict; kwargs...)
    end
    empty!(operator_dict)
    return new_operator_dict, operator_dict, param_idx
end


@inline function applygatetoone!(gate::PauliGateUnion, operator, coefficient, theta, param_idx, operator_dict, new_operator_dict, args...; kwargs...)

    if commutes(gate, operator)
        return operator_dict, new_operator_dict
    end

    coeff1 = applycos(coefficient, theta; param_idx=param_idx)
    sign, new_oper = getnewoperator(gate, operator)
    coeff2 = applysin(coefficient, theta; sign=sign, param_idx=param_idx)
    operator_dict[operator] = coeff1
    new_operator_dict[new_oper] = coeff2

    return
end

@inline function applygatetoone!(gate::CliffordGate, operator, coefficient, theta, param_idx, map_array, operator_dict, new_operator_dict, args...; kwargs...)
    # TODO: use more of the apply function here. This requires rewriting those.
    operator = copy(operator)
    qinds = gate.qinds
    lookup_op = _extractlookupop(operator, qinds)
    sign, new_op = map_array[lookup_op+1]  # +1 because Julia is 1-indexed and lookup_op is 0-indexed
    operator = _insertnewop!(operator, new_op, qinds)

    coefficient = _multiplysign!(coefficient, sign)

    new_operator_dict[operator] = coefficient

    return
end

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

function checktruncationonall!(operator_dict; max_weight::Real=Inf, min_abs_coeff=0.0, max_freq::Real=Inf, max_sins::Real=Inf, kwargs...)
    # TODO: This does currently hinder performance, even if we don't truncated
    # approx 55ms -> 66ms
    for (operator, coeff) in operator_dict
        checktruncationonone!(operator_dict, operator, coeff; max_weight=max_weight, min_abs_coeff=min_abs_coeff, max_freq=max_freq, max_sins=max_sins, kwargs...)
    end
    return
end

@inline function checktruncationonone!(operator_dict, operator, coeff; max_weight::Real=Inf, min_abs_coeff=0.0, max_freq::Real=Inf, max_sins::Real=Inf, kwargs...)
    we_truncate = false
    if truncateweight(operator, max_weight)
        we_truncate = true
    elseif truncatemincoeff(coeff, min_abs_coeff)
        we_truncate = true
    elseif truncatefrequency(coeff, max_freq)
        we_truncate = true
    elseif truncatesins(coeff, max_sins)
        we_truncate = true
    end
    if we_truncate
        delete!(operator_dict, operator)
    end
    return
end