### NOTE ###
# This is an experimental submodule of PauliPropagation.julia
# The Pauli propagation Surrogate works in the limited case of circuits consisting of Clifford gates and Pauli gates.
# The handling is awkward and not final.
# View this as legacy code.
# This code will be removed in the future.
############ 

"""
Abstract node type for the Pauli propagation Surrogate 
"""
abstract type CircuitNode end

"""
    EvalEndNode(operator::Integer, coefficient::Real)

Node type for the Pauli strings in the observable to be backpropagated.
"""
@kwdef mutable struct EvalEndNode <: CircuitNode
    operator::Integer
    coefficient::Real
    cummulative_value::Float64 = 0.0
    is_evaluated::Bool = false
end

"""
    EvalEndNode(operator::Integer)

Constructor for `EvalEndNode` with a default coefficient of 1.0.
"""
EvalEndNode(operator::Integer) = EvalEndNode(operator, 1.0)

"""
    PauliGateNode(parents::Vector{Union{EvalEndNode,PauliGateNode}}, trig_inds::Vector{Int}, signs::Vector{Int}, param_idx::Int)

Surrogate graph node for a Pauli gate.
"""
@kwdef mutable struct PauliGateNode <: CircuitNode
    parents::Vector{Union{EvalEndNode,PauliGateNode}}
    trig_inds::Vector{Int}
    signs::Vector{Int}
    param_idx::Int
    cummulative_value::Float64 = 0.0  # This must be changed to Real for automatic differentiation libraries.
    is_evaluated::Bool = false
end

"""
    NodePathProperties(coeff::CircuitNode, nsins::Int, ncos::Int, freq::Int)

Surrogate `PathProperties` type. Carries `CircuitNode`s instead of numerical coefficients.
"""
mutable struct NodePathProperties <: PathProperties
    coeff::CircuitNode
    nsins::Int
    ncos::Int
    freq::Int
end

"""
    NodePathProperties(coeff::CircuitNode)

One-argument constructor for `NodePathProperties`. Initializes `nsins`, `ncos`, and `freq` to 0.
"""
NodePathProperties(coeff::CircuitNode) = NodePathProperties(coeff, 0, 0, 0)

"""
    wrapcoefficients(pstr::PauliString, ::Type{NodePathProperties})

Wrap the coefficient of a `PauliString` into `NodePathProperties`.
"""
function wrapcoefficients(pstr::PauliString, ::Type{NodePathProperties})
    node = NodePathProperties(EvalEndNode(pstr.operator, pstr.coeff, 0.0, false))
    return PauliString(pstr.nqubits, pstr.operator, node)
end

"""
    wrapcoefficients(psum::PauliSum, ::Type{NodePathProperties})

Wrap the coefficients of a `PauliSum` into `NodePathProperties`.
"""
function wrapcoefficients(psum::PauliSum, ::Type{NodePathProperties})
    # node = NodePathProperties(EvalEndNode(pstr.operator, pstr.coeff, 0.0, false))
    # return PauliString(pstr.nqubits, pstr.operator, node)
    return PauliSum(psum.nqubits, Dict(op => NodePathProperties(EvalEndNode(op, coeff, 0.0, false)) for (op, coeff) in psum.op_dict))
end


function _multiplysign!(path_property::NodePathProperties, sign)
    _multiplysign!(path_property.coeff, sign)
    return path_property
end

function _multiplysign!(pauli_node::PauliGateNode, sign)
    for ii in eachindex(pauli_node.signs)
        pauli_node.signs[ii] *= sign
    end
    return pauli_node
end

function _multiplysign!(eval_endnode::EvalEndNode, sign)
    eval_endnode.coefficient *= sign
    return eval_endnode
end

## old
# function operatortopathdict(op, coefficient=1.0)
#     op = deepcopy(op)
#     d = Dict(op => NodePathProperties(EvalEndNode(operator=op, coefficient=coefficient, cummulative_value=coefficient)))
#     return d
# end

parents(node::T) where {T<:CircuitNode} = node.parents

function getnodeval(node::T) where {T<:CircuitNode}
    return node.cummulative_value
end

function isevaluated(node::T)::Bool where {T<:CircuitNode}
    return node.is_evaluated
end

function setcummulativevalue(node::CircuitNode, val)
    node.cummulative_value = val
    return
end

"""
Pretty print for `CircuitNode`
"""
Base.show(io::IO, node::CircuitNode) = print(io, "$(typeof(node))($(length(node.parents)) parent(s), param_idx=$(node.param_idx))")
"""
Pretty print for `EvalEndNode`
"""
Base.show(io::IO, node::EvalEndNode) = print(io, "$(typeof(node))(Operator=$(node.operator), coefficient=$(node.coefficient))")



function applycos(node::CircuitNode, theta; sign=1, param_idx=0)
    return PauliGateNode(parents=[node], trig_inds=[1], signs=[sign], param_idx=param_idx)
end

function applysin(node::CircuitNode, theta; sign=1, param_idx=0)
    return PauliGateNode(parents=[node], trig_inds=[-1], signs=[sign], param_idx=param_idx)
end

function merge(node1::CircuitNode, node2::CircuitNode)
    append!(node1.parents, node2.parents)
    append!(node1.trig_inds, node2.trig_inds)
    append!(node1.signs, node2.signs)
    return node1
end


function resetnodes(circuit_node::CircuitNode)
    if circuit_node.is_evaluated
        circuit_node.is_evaluated = false
        setcummulativevalue(circuit_node, 0.0)

        for parent_node in circuit_node.parents
            resetnodes(parent_node)
        end
    end
end

function resetnodes(end_node::EvalEndNode)
    end_node.is_evaluated = false
    setcummulativevalue(end_node, 0.0)
    return
end

function resetnodes(eval_list::Vector{<:CircuitNode})
    @threads for node in eval_list
        node.is_evaluated = false
    end
    return
end

"""
    expectation(eval_list::Vector{<:CircuitNode}, thetas)

Evaluate the expectation value of a Surrogate by evaluating all involved circuit nodes in the correct order.
`eval_list` can be attained as the output of `gettraceevalorder()`
"""
function expectation(eval_list::Vector{<:CircuitNode}, thetas)
    for ii in eachindex(eval_list)
        _evalnode(eval_list[ii], thetas)
    end
    return getnodeval(eval_list[end])
end

"""
    gettraceevalorder(node::CircuitNode, thetas)

Return a vector of `CircuitNode`s in the order they should be evaluated to get the correct cummulative result on `node`.
`thetas` numerically plays no role here but it needs to be the correct length given the number of parametrized gates.
"""
function gettraceevalorder(node::CircuitNode, thetas)
    eval_list = Union{PauliGateNode,EvalEndNode}[]
    traceevalorder(node, thetas; eval_list)
    return eval_list
end

"""
    traceevalorder(node::PauliGateNode, thetas; eval_list=nothing)

Evaluate the coefficient of `node` on the Surrogate by recursively evaluating all parents.
`thetas` are the parameters of the circuit.
NOTE: This requires calling `resetnodes` in-between evaluations with different `thetas`.
`eval_list` does not need to be passed when manually using this function.
"""
function traceevalorder(node::PauliGateNode, thetas; eval_list=nothing)
    val = typeof(node.cummulative_value)(0)

    if isevaluated(node)
        return node.cummulative_value
    end

    node_parents = parents(node)

    for ii in eachindex(node_parents)
        parent = node_parents[ii]

        this_val = _evaltrig(node.trig_inds[ii], node.signs[ii], thetas, node.param_idx)

        other_val = isevaluated(parent) ? getnodeval(parent) : traceevalorder(parent, thetas; eval_list)
        val += this_val * other_val

    end
    setcummulativevalue(node, val)
    node.is_evaluated = true
    if !isnothing(eval_list)
        push!(eval_list, node)
    end
    return node.cummulative_value
end

"""
    traceevalorder(node::EvalEndNode, thetas; eval_list=nothing)

Evaluates the observable's coefficient. 
This function likely does not need to be called manually.
"""
function traceevalorder(node::EvalEndNode, thetas; eval_list=nothing)
    setcummulativevalue(node, node.coefficient)
    node.is_evaluated = true
    if !isnothing(eval_list)
        push!(eval_list, node)
    end
    return node.cummulative_value
end

"""
    traceevalorder(nodes::Vector{<:CircuitNode}, thetas)

Evaluate the sum of coefficients of a vector of `CircuitNode` on the Surrogate. 
This will be evaluated in parallel with recursive evaluation of the parents. 
NOTE: This requires calling `resetnodes` in-between evaluations with different `thetas`.
"""
traceevalorder(nodes::Vector{<:CircuitNode}, thetas) = ThreadsX.sum(traceevalorder(node, thetas) for node in nodes)


# const trig_funcs = [x -> 1.0, cos, sin] 

function _evaltrig(which_idx, sign, thetas, param_idx)
    if which_idx == 1
        return cos(thetas[param_idx]) * sign
    elseif which_idx == -1
        return sin(thetas[param_idx]) * sign
    else
        return 1.0 * sign
    end
    # return trig_funcs[which_idx+1](thetas[param_idx]) * sign
end

_evaltrig(node::CircuitNode, thetas, ii) = _evaltrig(node.trig_inds[ii], node.signs[ii], thetas, node.param_idx)



## Exploratory

function _evalnode(node::PauliGateNode, thetas)

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


# TODO: not functional
# """
#     toevallist(final_nodes::Vector{<:CircuitNode})

# Return a vector of `CircuitNode`s in the order they should be evaluated to get the correct cummulative result on the final node.
# """
# function toevallist(final_nodes::Vector{<:CircuitNode})
#     final_eval_node = PauliGateNode(
#         parents=final_nodes,
#         trig_inds=zeros(Int, length(final_nodes)),
#         signs=ones(length(final_nodes)),
#         param_idx=1,
#         cummulative_value=0.0
#     )

#     resetnodes(final_eval_node)
#     resetnodes(final_eval_node)
#     return gettraceevalorder(final_eval_node, zeros(m_total))
# end