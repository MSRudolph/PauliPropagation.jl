### Propagation necessities

"""
    propagate(circ, pstr::PauliString{PauliStringType,NodePathProperties}; max_weight=Inf, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, kwargs...)

Construct a Pauli propagation surrogate of the propagated `PauliString` through the circuit `circ` in the Heisenberg picture. 
The circuit must only contain `CliffordGate`s and `PauliRotation`s.
It is applied to the Pauli string in reverse order, and the action of each gate is its conjugate action.
Default truncations are `max_weight`, `max_freq`, and `max_sins`.
A custom truncation function can be passed as `customtruncfunc` with the signature customtruncfunc(pstr::PauliStringType, coefficient)::Bool.
Further `kwargs` are passed to the lower-level functions `applymergetruncate!`, `applytoall!`, `applyandadd!`, and `apply`.
"""
function propagate(circ, pstr::PauliString{TT,NodePathProperties}; max_weight=Inf, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, kwargs...) where {TT<:PauliStringType}
    return propagate(circ, PauliSum(pstr); max_weight, max_freq, max_sins, customtruncfunc, kwargs...)
end

"""
    propagate(circ, psum::PauliSum{PauliStringType,NodePathProperties}; max_weight=Inf, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, kwargs...)

Construct a Pauli propagation surrogate of the propagated `PauliSum` through the circuit `circ` in the Heisenberg picture.
The circuit must only contain `CliffordGate`s and `PauliRotation`s. 
It is applied to the Pauli sum in reverse order, and the action of each gate is its conjugate action.
Default truncations are `max_weight`, `max_freq`, and `max_sins`.
A custom truncation function can be passed as `customtruncfunc` with the signature customtruncfunc(pstr::PauliStringType, coefficient)::Bool.
Further `kwargs` are passed to the lower-level functions `applymergetruncate!`, `applytoall!`, `applyandadd!`, and `apply`.
"""
function propagate(circ, psum::PauliSum{TT,NodePathProperties}; max_weight=Inf, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, kwargs...) where {TT<:PauliStringType}
    _checksurrogationconditions(circ)
    return propagate!(circ, PauliSum(psum.nqubits, copy(psum.terms)); max_weight, max_freq, max_sins, customtruncfunc, kwargs...)
end

"""
    propagate!(circ, psum::PauliSum{PauliStringType,NodePathProperties}; max_weight=Inf, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, kwargs...)

Construct a Pauli propagation surrogate of the propagated `PauliSum` through the circuit `circ` in the Heisenberg picture. 
The `PauliSum` `psum` is modified in place.
The circuit must only contain `CliffordGate`s and `PauliRotation`s. 
It is applied to the Pauli sum in reverse order, and the action of each gate is its conjugate action.
The input `psum` will be modified.
Default truncations are `max_weight`, `max_freq`, and `max_sins`.
A custom truncation function can be passed as `customtruncfunc` with the signature customtruncfunc(pstr::PauliStringType, coefficient)::Bool.
Further `kwargs` are passed to the lower-level functions `applymergetruncate!`, `applytoall!`, `applyandadd!`, and `apply`.
"""
function propagate!(circ, psum::PauliSum{TT,NodePathProperties}; max_weight=Inf, max_freq=Inf, max_sins=Inf, customtruncfunc=nothing, kwargs...) where {TT<:PauliStringType}
    _checksurrogationconditions(circ)

    thetas = Array{Float64}(undef, countparameters(circ))

    param_idx = length(thetas)

    aux_psum = similar(psum)

    for gate in reverse(circ)
        # add param_index as kwarg, which will descend into the apply function eventually
        psum, aux_psum, param_idx = applymergetruncate!(gate, psum, aux_psum, thetas, param_idx; max_weight, max_freq, max_sins, customtruncfunc, param_idx=param_idx, kwargs...)
    end
    return psum
end

function _checksurrogationconditions(circ)
    if !all(isa(gate, CliffordGate) || isa(gate, PauliRotation) for gate in circ)
        throw(ArgumentError("The surrogate currently only accepts `CliffordGate`s and `PauliRotation`s."))
    end
    return
end


## For Pauli Rotations
"""
    applytoall!(gate::PauliRotation, theta, psum, aux_psum; kwargs...)

Overload of `applytoall!` for `PauliRotation` gates. 
It fixes the type-instability of the `apply()` function and reduces moving Pauli strings between `psum` and `aux_psum`.
`psum` and `aux_psum` are merged later.
"""
function applytoall!(gate::PauliRotation, theta, psum::PauliSum{TT,NodePathProperties}, aux_psum; kwargs...) where {TT<:PauliStringType}
    # turn the (potentially) PauliRotation gate into a MaskedPauliRotation gate
    # this allows for faster operations
    gate = _tomaskedpaulirotation(gate, paulitype(psum))

    # loop over all Pauli strings and their coefficients in the Pauli sum
    for (pstr, coeff) in psum

        if commutes(gate, pstr)
            # if the gate commutes with the pauli string, do nothing
            continue
        end

        # else we know the gate will split th Pauli string into two
        pstr, coeff1, new_pstr, coeff2 = splitapply(gate, pstr, coeff, theta; kwargs...)

        # set the coefficient of the original Pauli string
        set!(psum, pstr, coeff1)

        # set the coefficient of the new Pauli string in the aux_psum
        # we can set the coefficient because PauliRotations create non-overlapping new Pauli strings
        set!(aux_psum, new_pstr, coeff2)
    end

    return
end

function splitapply(gate::MaskedPauliRotation, pstr::PauliStringType, coeff::NodePathProperties, theta; kwargs...)
    coeff1 = _applycos(coeff, theta; kwargs...)
    new_pstr, sign = getnewpaulistring(gate, pstr)
    coeff2 = _applysin(coeff, theta, sign; kwargs...)

    return pstr, coeff1, new_pstr, coeff2
end

function _applycos(path::NodePathProperties, theta, sign=1; param_idx=0, kwargs...)
    return NodePathProperties(_applycos(path.node, theta, sign; param_idx=param_idx), path.nsins, path.ncos + 1, path.freq + 1)
end

function _applycos(node::CircuitNode, theta, sign=1; param_idx=0, kwargs...)
    return PauliRotationNode(parents=[node], trig_inds=[1], signs=[sign], param_idx=param_idx)
end

function _applysin(path::NodePathProperties, theta, sign=1; param_idx=0, kwargs...)
    return NodePathProperties(_applysin(path.node, theta, sign; param_idx=param_idx), path.nsins + 1, path.ncos, path.freq + 1)
end

function _applysin(node::CircuitNode, theta, sign=1; param_idx=0, kwargs...)
    return PauliRotationNode(parents=[node], trig_inds=[-1], signs=[sign], param_idx=param_idx)
end

function merge(pth1::NodePathProperties, pth2::NodePathProperties)
    return NodePathProperties(
        merge(pth1.node, pth2.node),
        min(pth1.nsins, pth2.nsins),
        min(pth1.ncos, pth2.ncos),
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

"""
    apply(gate::CliffordGate, pstr::PauliStringType, coeff::NodePathProperties)

Apply a `CliffordGate` to an integer Pauli string and `NodePathProperties` coefficient. 
"""
function apply(gate::CliffordGate, pstr::PauliStringType, coeff::NodePathProperties; kwargs...)
    # this array carries the new Paulis + sign for every occuring old Pauli combination
    map_array = clifford_map[gate.symbol]

    qinds = gate.qinds

    # this integer carries the active Paulis on its bits
    lookup_int = getpauli(pstr, qinds)

    # this integer can be used to index into the array returning the new Paulis
    # +1 because Julia is 1-indexed and lookup_int is 0-indexed
    sign, partial_pstr = map_array[lookup_int+1]

    # insert the bits of the new Pauli into the old Pauli
    pstr = setpauli(pstr, partial_pstr, qinds)

    coeff = _multiplysign(coeff, sign)

    return pstr, coeff
end

function _multiplysign(pth::NodePathProperties, sign; kwargs...)
    return NodePathProperties(_multiplysign(pth.node, sign), pth.nsins, pth.ncos, pth.freq)
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