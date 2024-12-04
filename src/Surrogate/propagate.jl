### Propagation necessities

"""
    propagate(circ, pstr::PauliString; kwargs...)

Propagate a `PauliString` through the circuit `circ` in the Heisenberg picture. 
`kwargs` are passed to the truncation function. Supported by default are `max_weight`, `min_abs_coeff`, `max_freq`, and `max_sins`.
A custom truncation function can be passed as `customtruncatefn` with the signature customtruncatefn(pstr::PauliStringType, coefficient)::Bool.
"""
function propagate(circ, pstr::PauliString; kwargs...)
    _checkcoefftype(pstr)
    psum = PauliSum(pstr.nqubits, pstr)
    return propagate(circ, psum; kwargs...)
end

"""
    propagate(circ, psum::PauliSum; kwargs...)

Propagate a `PauliSum` through the circuit `circ` in the Heisenberg picture. 
`kwargs` are passed to the truncation function. Supported by default are `max_weight`, `min_abs_coeff`, `max_freq`, and `max_sins`.
A custom truncation function can be passed as `customtruncatefn` with the signature customtruncatefn(pstr::PauliStringType, coefficient)::Bool.
"""
function propagate(circ, psum::PauliSum; kwargs...)
    _checkcoefftype(psum)
    pauli_dict = propagate!(circ, deepcopy(psum.op_dict); kwargs...)
    return PauliSum(psum.nqubits, pauli_dict)
end


"""
    propagate!(circ, psum::PauliSum; kwargs...)

Propagate a `PauliSum` through the circuit `circ` in the Heisenberg picture. 
The input `psum` will be modified.
`kwargs` are passed to the truncation function. Supported by default are `max_weight`, `min_abs_coeff`, `max_freq`, and `max_sins`.
A custom truncation function can be passed as `customtruncatefn` with the signature customtruncatefn(pstr::PauliStringType, coefficient)::Bool.
"""
function propagate!(circ, psum::PauliSum; kwargs...)
    # check that circ only constists of Pauli gates and Clifford gates
    if !all(isa(gate, CliffordGate) || isa(gate, PauliGateUnion) for gate in circ)
        throw(ArgumentError("The surrogate currently only accepts Clifford gates and (Fast)Pauli gates."))
    end

    propagate!(circ, psum.op_dict; kwargs...)
    return psum
end

"""
    propagate!(circ, d::Dict{TermType,NodePathProperties}; kwargs...)

Propagate a `Dict{PauliStringType,CoeffType}` through the circuit `circ` in the Heisenberg picture. 
The input `psum` will be modified.
`kwargs` are passed to the truncation function. Supported by default are `max_weight`, `min_abs_coeff`, `max_freq`, and `max_sins`.
A custom truncation function can be passed as `customtruncatefn` with the signature customtruncatefn(pstr::PauliStringType, coefficient)::Bool.
"""
function propagate!(circ, psum::Dict{TermType,NodePathProperties}; kwargs...) where {TermType<:PauliStringType}
    thetas = Array{Float64}(undef, countparameters(circ))

    param_idx = length(thetas)

    second_psum = typeof(psum)()  # pre-allocating somehow doesn't do anything

    ## TODO:
    # - decide where to reverse the circuit
    # - verbose option  
    # - more elegant param_idx incrementation
    for gate in reverse(circ)
        # add param_index as kwarg, which will descend into the apply function eventually
        psum, second_psum, param_idx = mergingapply!(gate, psum, second_psum, thetas, param_idx; param_idx=param_idx, kwargs...)
    end
    return psum
end


function _checkcoefftype(pobj::Union{PauliString,PauliSum})
    # convert numerical coefficients to `NodePathProperties` 
    CoeffType = typeof(pobj).parameters[2]
    if CoeffType <: Number
        throw(
            "`You are using the Surrogate's propagation function without passing parameters. " *
            "But the current coefficient type is $(CoeffType)`. " *
            "Consider converting to `NodePathProperties` via `wrapcoefficients(your_current_paulis, NodePathProperties)`."
        )
    end
    return
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