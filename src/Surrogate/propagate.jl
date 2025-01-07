### Propagation necessities

"""
    propagate(circ, pstr::PauliString{PauliStringType,NodePathProperties}; kwargs...)

Construct a Pauli propagation surrogate of the propagated `PauliString` through the circuit `circ` in the Heisenberg picture. 
The circuit must only contain `CliffordGate`s and `PauliRotation`s.
It is applied to the Pauli string in reverse order, and the action of each gate is its conjugate action.
`kwargs` are passed to the truncation function. Supported by default for the surrogation are `max_weight`, `max_freq`, and `max_sins`.
A custom truncation function can be passed as `customtruncatefn` with the signature customtruncatefn(pstr::PauliStringType, coefficient)::Bool.
"""
function propagate(circ, pstr::PauliString{TT,NodePathProperties}; kwargs...) where {TT<:PauliStringType}
    return propagate(circ, PauliSum(pstr); kwargs...)
end

"""
    propagate(circ, psum::PauliSum{PauliStringType,NodePathProperties}; kwargs...)

Construct a Pauli propagation surrogate of the propagated `PauliSum` through the circuit `circ` in the Heisenberg picture.
The circuit must only contain `CliffordGate`s and `PauliRotation`s. 
It is applied to the Pauli sum in reverse order, and the action of each gate is its conjugate action.
`kwargs` are passed to the truncation function. Supported by default for the surrogation are `max_weight`, `max_freq`, and `max_sins`.
A custom truncation function can be passed as `customtruncatefn` with the signature customtruncatefn(pstr::PauliStringType, coefficient)::Bool.
"""
function propagate(circ, psum::PauliSum{TT,NodePathProperties}; kwargs...) where {TT<:PauliStringType}
    _checksurrogationconditions(circ)
    return propagate!(circ, PauliSum(psum.nqubits, copy(psum.terms)); kwargs...)
end

"""
    propagate!(circ, psum::PauliSum{PauliStringType,NodePathProperties}; kwargs...)

Construct a Pauli propagation surrogate of the propagated `PauliSum` through the circuit `circ` in the Heisenberg picture. 
The `PauliSum` `psum` is modified in place.
The circuit must only contain `CliffordGate`s and `PauliRotation`s. 
It is applied to the Pauli sum in reverse order, and the action of each gate is its conjugate action.
The input `psum` will be modified.
`kwargs` are passed to the truncation function. Supported by default for the surrogation are `max_weight`, `max_freq`, and `max_sins`.
A custom truncation function can be passed as `customtruncatefn` with the signature customtruncatefn(pstr::PauliStringType, coefficient)::Bool.
"""
function propagate!(circ, psum::PauliSum{TT,NodePathProperties}; kwargs...) where {TT<:PauliStringType}
    _checksurrogationconditions(circ)

    thetas = Array{Float64}(undef, countparameters(circ))

    param_idx = length(thetas)

    aux_psum = similar(psum)

    for gate in reverse(circ)
        # add param_index as kwarg, which will descend into the apply function eventually
        psum, aux_psum, param_idx = applymergetruncate!(gate, psum, aux_psum, thetas, param_idx; param_idx=param_idx, kwargs...)
    end
    return psum
end

function _checksurrogationconditions(circ)
    if !all(isa(gate, CliffordGate) || isa(gate, PauliRotation) for gate in circ)
        throw(ArgumentError("The surrogate currently only accepts `CliffordGate`s and `PauliRotation`s."))
    end
    return
end


## For Pauli Gates

function splitapply(gate::MaskedPauliRotation, pstr::PauliStringType, coeff::NodePathProperties, theta; kwargs...)
    coeff1 = _applycos(coeff, theta; kwargs...)
    new_pstr, sign = getnewpaulistring(gate, pstr)
    coeff2 = _applysin(coeff, theta, sign; kwargs...)

    return pstr, coeff1, new_pstr, coeff2
end

function _applycos(node::CircuitNode, theta, sign=1; param_idx=0, kwargs...)
    return PauliRotationNode(parents=[node], trig_inds=[1], signs=[sign], param_idx=param_idx)
end

function _applysin(node::CircuitNode, theta, sign=1; param_idx=0, kwargs...)
    return PauliRotationNode(parents=[node], trig_inds=[-1], signs=[sign], param_idx=param_idx)
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

function _multiplysign(pauli_node::PauliRotationNode, sign; kwargs...)
    for ii in eachindex(pauli_node.signs)
        pauli_node.signs[ii] *= sign
    end
    return pauli_node
end

function _multiplysign(eval_endnode::EvalEndNode, sign; kwargs...)
    eval_endnode.coefficient *= sign
    return eval_endnode
end

## Truncation functions

# don't truncate on coefficients
function truncatemincoeff(path::NodePathProperties, min_abs_coeff::Real)
    return false
end