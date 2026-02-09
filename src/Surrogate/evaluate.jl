"""
    evaluate!(psum::PauliSum{<:Integer,NodePathProperties}, thetas; reset=true)

Evaluate the expectation value of a Surrogate by evaluating all involved circuit nodes in the correct order.
`eval_list` can be attained as the output of `gettraceevalorder()`.
If `reset` is false, the function will not reset the `is_evaluated` flags of the nodes.
This should only be done if they are manually reset.
"""
function evaluate!(psum::PauliSum{TT,NodePathProperties}, thetas; reset=true) where {TT}
    # collect the final paths first for better multi-threaded iteration
    paths = collect(coefficients(psum))

    # set all is_evaluated flags to false
    if reset
        reset!(paths)
    end

    # eval recursively
    @threads for pth in paths
        PauliPropagation._traceevalorder(pth.node, thetas)
    end

    return psum
end


## Reset functions
"""
    reset!(psum::PauliSum{<:Integer, NodePathProperties})

Reset the nodes in a the Surrogate. 
Needs to be done in-between evaluations with different parameters.
"""
function reset!(psum::PauliSum{TT,CT}) where {TT,CT<:NodePathProperties}
    paths = collect(coefficients(psum))
    reset!(paths)
    return
end

"""
    reset!(paths::Vector{NodePathProperties}) 

Reset a vector of `NodePathProperties` in a the Surrogate. 
Needs to be done in-between evaluations with different parameters.
"""
function reset!(paths::AbstractVector{CT}) where {CT<:NodePathProperties}
    @threads for pth in paths
        reset!(pth.node)
    end
    return
end

"""
    reset!(circuit_node::CircuitNode)

Reset a `CircuitNode` in a the Surrogate. Needs to be done in-between evaluations with different parameters.
"""
function reset!(circuit_node::CircuitNode)
    if circuit_node.is_evaluated
        circuit_node.is_evaluated = false
        setcummulativevalue(circuit_node, 0.0)

        for parent_node in circuit_node.parents
            reset!(parent_node)
        end
    end
end

"""
    reset!(end_node::EvalEndNode)

Reset a `EvalEndNode` in a the Surrogate. 
These sit on the `coeff` field of `NodePathProperties` and are the end of the evaluation chain.
Needs to be done in-between evaluations with different parameters.
"""
function reset!(end_node::EvalEndNode)
    end_node.is_evaluated = false
    setcummulativevalue(end_node, 0.0)
    return
end

## Evaluation functions

# Evaluate the coefficient of `node` on the Surrogate by recursively evaluating all parents.
# `thetas` are the parameters of the circuit.
# NOTE: This requires calling `resetnodes` in-between evaluations with different `thetas`.
# `eval_list` does not need to be passed when manually using this function.
function _traceevalorder(node::PauliRotationNode, thetas; eval_list=nothing)
    val = zero(node.cummulative_value)

    if isevaluated(node)
        return node.cummulative_value
    end

    node_parents = parents(node)

    for ii in eachindex(node_parents)
        parent = node_parents[ii]

        this_val = _evaltrig(node.trig_inds[ii], node.signs[ii], thetas, node.param_idx)

        other_val = isevaluated(parent) ? getnodeval(parent) : _traceevalorder(parent, thetas; eval_list)
        val += this_val * other_val

    end
    setcummulativevalue(node, val)
    node.is_evaluated = true
    if !isnothing(eval_list)
        push!(eval_list, node)
    end
    return node.cummulative_value
end


# Evaluates the observable's coefficient. 
# This function likely does not need to be called manually.
function _traceevalorder(node::EvalEndNode, thetas; eval_list=nothing)
    setcummulativevalue(node, node.coefficient)
    node.is_evaluated = true
    if !isnothing(eval_list)
        push!(eval_list, node)
    end
    return node.cummulative_value
end


# Evaluate the sum of coefficients of a vector of `CircuitNode` on the Surrogate. 
# This will be evaluated in parallel with recursive evaluation of the parents. 
# NOTE: This requires calling `resetnodes` in-between evaluations with different `thetas`.
_traceevalorder(nodes::Vector{<:CircuitNode}, thetas) = sum(_traceevalorder(node, thetas) for node in nodes)


# when evaluating the surrogate, we have values stored on the edges indicating a sine or cosine coefficient
function _evaltrig(which_idx, sign, thetas, param_idx)
    if which_idx == 1
        return cos(thetas[param_idx]) * sign
    elseif which_idx == -1
        return sin(thetas[param_idx]) * sign
    else
        return 1.0 * sign
    end
end

_evaltrig(node::CircuitNode, thetas, ii) = _evaltrig(node.trig_inds[ii], node.signs[ii], thetas, node.param_idx)


function _evalnode(node::PauliRotationNode, thetas)

    val = sum(_evaltrig(node, thetas, ii) * getnodeval(node.parents[ii]) for ii in eachindex(node.parents))
    # val = sum(evaltrig(node.trig_inds[ii], node.signs[ii], numeric_thetas, node.param_idx) * getnodeval(node.parents[ii]) for ii in eachindex(node.parents))
    setcummulativevalue(node, val)
    node.is_evaluated = true
    return
end

function _evalnode(node::EvalEndNode, thetas)
    setcummulativevalue(node, node.coefficient)
    node.is_evaluated = true
    return
end



### Area for connoisseurs


# Evaluate the expectation value of a Surrogate by evaluating all involved circuit nodes in the correct order.
# `eval_list` can be attained as the output of `gettraceevalorder()`
function expectation(eval_list::Vector{<:CircuitNode}, thetas)
    for ii in eachindex(eval_list)
        _evalnode(eval_list[ii], thetas)
    end
    return getnodeval(eval_list[end])
end


# Return a vector of `CircuitNode`s in the order they should be evaluated to get the correct cummulative result on `node`.
# `thetas` numerically plays no role here but it needs to be the correct length given the number of parametrized gates.
function gettraceevalorder(node::CircuitNode, thetas)
    eval_list = Union{PauliRotationNode,EvalEndNode}[]
    _traceevalorder(node, thetas; eval_list)
    return eval_list
end


# Resets a vector of `CircuitNode`s.
function reset!(eval_list::Vector{<:CircuitNode})
    @threads for node in eval_list
        node.is_evaluated = false
    end
    return
end