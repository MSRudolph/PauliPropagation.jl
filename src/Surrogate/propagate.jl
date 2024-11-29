### Mergingbfs necessities

"""
    mergingbfs(circ, pstr::PauliString; kwargs...)

Perform merging breadth-first search (BFS) surrogation of a `PauliString` propagating through the circuit `circ` in the Heisenberg picture. 
`kwargs` are passed to the truncation function. Supported by default are `max_weight`, `min_abs_coeff`, `max_freq`, and `max_sins`.
A custom truncation function can be passed as `customtruncatefn` with the signature customtruncatefn(pstr::PauliStringType, coefficient)::Bool.
"""
function mergingbfs(circ, pstr::PauliString; kwargs...)
    pstr = _checkandconvert(pstr)
    psum = PauliSum(pstr.nqubits, pstr)
    return mergingbfs(circ, psum; kwargs...)
end

"""
    mergingbfs(circ, psum::PauliSum; kwargs...)

Perform merging breadth-first search (BFS) surrogation of a `PauliSum` propagating through the circuit `circ` in the Heisenberg picture. 
`kwargs` are passed to the truncation function. Supported by default are `max_weight`, `min_abs_coeff`, `max_freq`, and `max_sins`.
A custom truncation function can be passed as `customtruncatefn` with the signature customtruncatefn(pstr::PauliStringType, coefficient)::Bool.
"""
function mergingbfs(circ, psum::PauliSum; kwargs...)
    psum = _checkandconvert(psum)
    pauli_dict = mergingbfs!(circ, deepcopy(psum.op_dict); kwargs...)
    return PauliSum(psum.nqubits, pauli_dict)
end


"""
    mergingbfs!(circ, psum::PauliSum; kwargs...)

Perform in-place merging breadth-first search (BFS) surrogation of a `PauliSum` propagating through the circuit `circ` in the Heisenberg picture. 
The input `psum` will be modified.
`kwargs` are passed to the truncation function. Supported by default are `max_weight`, `min_abs_coeff`, `max_freq`, and `max_sins`.
A custom truncation function can be passed as `customtruncatefn` with the signature customtruncatefn(pstr::PauliStringType, coefficient)::Bool.
"""
function mergingbfs!(circ, psum::PauliSum; kwargs...)
    # check that circ only constists of Pauli gates and Clifford gates
    if !all(isa(gate, CliffordGate) || isa(gate, PauliGateUnion) for gate in circ)
        throw(ArgumentError("The surrogate currently only accepts Clifford gates and (Fast)Pauli gates."))
    end

    pauli_dict = mergingbfs!(circ, psum.op_dict; kwargs...)
    return PauliSum(psum.nqubits, pauli_dict)
end

"""
    mergingbfs!(circ, d::Dict; kwargs...)

Perform in-place merging breadth-first search (BFS) simulation of a `Dict{PauliStringType,CoeffType}` propagating through the circuit `circ` in the Heisenberg picture. 
The input `psum` will be modified.
`kwargs` are passed to the truncation function. Supported by default are `max_weight`, `min_abs_coeff`, `max_freq`, and `max_sins`.
A custom truncation function can be passed as `customtruncatefn` with the signature customtruncatefn(pstr::PauliStringType, coefficient)::Bool.
"""
function mergingbfs!(circ, d::Dict{OpType,NodePathProperties}; kwargs...) where {OpType<:PauliStringType}
    thetas = Array{Float64}(undef, countparameters(circ))
    # TODO: Should we have a version of this that doesn't require thetas and uses surrogate code?
    param_idx = length(thetas)

    second_d = typeof(d)()  # pre-allocating somehow doesn't do anything

    ## TODO:
    # - decide where to reverse the circuit
    # - verbose option  
    # - more elegant param_idx incrementation
    for gate in reverse(circ)
        # add param_index as kwarg, which will descend into the apply function eventually
        d, second_d, param_idx = mergingapply(gate, d, second_d, thetas, param_idx; param_idx=param_idx, kwargs...)
    end
    return d
end


function _checkandconvert(pobj::Union{PauliString,PauliSum})
    # convert numerical coefficients to `NodePathProperties` 
    CoeffType = typeof(pobj).parameters[2]
    if CoeffType <: Number
        println("Converting coefficient of type `$(CoeffType)` to `NodePathProperties` for the Surrogate.")
        pobj = wrapcoefficients(pobj, NodePathProperties)
    end
    return pobj
end


## For Pauli Gates

function applycos(node::CircuitNode, theta; sign=1, param_idx=0, kwargs...)
    return PauliGateNode(parents=[node], trig_inds=[1], signs=[sign], param_idx=param_idx)
end

function applysin(node::CircuitNode, theta; sign=1, param_idx=0, kwargs...)
    return PauliGateNode(parents=[node], trig_inds=[-1], signs=[sign], param_idx=param_idx)
end

function merge(pth1::NodePathProperties, pth2::NodePathProperties)
    return NodePathProperties(
        merge(pth1.coeff, pth2.coeff),
        min(pth1.ncos, pth2.ncos),
        min(pth1.nsins, pth2.nsins),
        min(pth1.freq, pth2.freq)
    )
end

function merge(node1::CircuitNode, node2::CircuitNode)
    append!(node1.parents, node2.parents)
    append!(node1.trig_inds, node2.trig_inds)
    append!(node1.signs, node2.signs)
    return node1
end

## For Clifford Gates

function _multiplysign(pth::NodePathProperties, sign; kwargs...)
    return NodePathProperties(_multiplysign(pth.coeff, sign), pth.nsins, pth.ncos, pth.freq)
end

function _multiplysign(pauli_node::PauliGateNode, sign; kwargs...)
    for ii in eachindex(pauli_node.signs)
        pauli_node.signs[ii] *= sign
    end
    return pauli_node
end

function _multiplysign(eval_endnode::EvalEndNode, sign; kwargs...)
    eval_endnode.coefficient *= sign
    return eval_endnode
end