## For the surrogate
abstract type CircuitNode end


@kwdef mutable struct EvalEndNode <: CircuitNode
    operator::Integer
    coefficient::Real # eventually symbolic?
    # func::Function = x->1.0
    cummulative_value::Float64 = 0.0
    is_evaluated::Bool = false

end

EvalEndNode(operator) = EvalEndNode(operator, 1.0)

@kwdef mutable struct PauliGateNode <: CircuitNode #where {T<:Real}
    parents::Union{Vector{EvalEndNode},Vector{PauliGateNode}}
    trig_inds::Vector{Int}
    signs::Vector{Int}
    param_idx::Int
    cummulative_value::Float64 = 0.0
    is_evaluated::Bool = false
end

mutable struct NodePathProperties <: PathProperties
    coeff::CircuitNode
    nsins::Int
    ncos::Int
    freq::Int
end

NodePathProperties(coeff) = NodePathProperties(coeff, 0, 0, 0)

function wrapcoefficients(pstr::PauliString, ::Type{NodePathProperties})
    node = NodePathProperties(EvalEndNode(pstr.operator, pstr.coeff, 0.0, false))
    return PauliString(pstr.nqubits, pstr.operator, node)
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

function operatortopathdict(op, coefficient=1.0)
    op = deepcopy(op)
    # d = Dict(op => PathProperties(EvalEndNode(operator=op, coefficient=coefficient, cummulative_value=coefficient)))
    d = Dict(op => NodePathProperties(EvalEndNode(operator=op, coefficient=coefficient, cummulative_value=coefficient)))
    return d
end

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

Base.show(io::IO, node::CircuitNode) = print(io, "$(typeof(node))($(length(node.parents)) parent(s), param_idx=$(node.param_idx))")
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


function gettraceevalorder(node::CircuitNode, numeric_thetas)
    eval_list = Union{PauliGateNode,EvalEndNode}[]
    traceevalorder(node, numeric_thetas; eval_list)
    return eval_list
end

function traceevalorder(node::PauliGateNode, numeric_thetas; eval_list=nothing)
    val = 0.0

    if isevaluated(node)
        return node.cummulative_value
    end

    node_parents = parents(node)

    for ii in eachindex(node_parents)
        # exprfunc = node.funcs[ii]
        parent = node_parents[ii]

        this_val = evaltrig(node.trig_inds[ii], node.signs[ii], numeric_thetas, node.param_idx)

        other_val = isevaluated(parent) ? getnodeval(parent) : traceevalorder(parent, numeric_thetas; eval_list)
        val += this_val * other_val

    end
    setcummulativevalue(node, val)
    node.is_evaluated = true
    if !isnothing(eval_list)
        push!(eval_list, node)
    end
    return node.cummulative_value
end

function traceevalorder(node::EvalEndNode, numeric_thetas; eval_list=nothing)
    setcummulativevalue(node, node.coefficient) #node.func(numeric_thetas)
    node.is_evaluated = true
    if !isnothing(eval_list)
        push!(eval_list, node)
    end
    return node.cummulative_value
end

traceevalorder(nodes::Vector{<:CircuitNode}, num_thetas) = ThreadsX.sum(traceevalorder(node, num_thetas) for node in nodes)


# const trig_funcs = [x -> 1.0, cos, sin] 

function evaltrig(which_idx, sign, thetas, param_idx)
    if which_idx == 1
        return cos(thetas[param_idx]) * sign
    elseif which_idx == -1
        return sin(thetas[param_idx]) * sign
    else
        return 1.0 * sign
    end
    # return trig_funcs[which_idx+1](thetas[param_idx]) * sign
end

evaltrig(node::CircuitNode, thetas, ii) = evaltrig(node.trig_inds[ii], node.signs[ii], thetas, node.param_idx)



## Exploratory

function evalnode(node::PauliGateNode, numeric_thetas)

    val = sum(evaltrig(node, numeric_thetas, ii) * getnodeval(node.parents[ii]) for ii in eachindex(node.parents))
    # val = sum(evaltrig(node.trig_inds[ii], node.signs[ii], numeric_thetas, node.param_idx) * getnodeval(node.parents[ii]) for ii in eachindex(node.parents))
    setcummulativevalue(node, val)
    node.is_evaluated = true
    return
end

function evalnode(node::EvalEndNode, numeric_thetas)
    setcummulativevalue(node, node.coefficient)
    node.is_evaluated = true
    return
end

function expectation(eval_list::Vector{<:CircuitNode}, numeric_thetas)
    for ii in eachindex(eval_list)
        evalnode(eval_list[ii], numeric_thetas)
    end
    return getnodeval(eval_list[end])
end


function toevallist(final_nodes)
    # final_nodes = collect(pth.coeff for (obs, pth) in dsym if !containsYorZ(obs));
    final_eval_node = PauliGateNode(parents=final_nodes, trig_inds=zeros(Int, length(final_nodes)), signs=ones(length(final_nodes)), param_idx=1, cummulative_value=0.0)

    resetnodes(final_eval_node)
    resetnodes(final_eval_node)
    return eval_list = gettraceevalorder(final_eval_node, zeros(m_total))
end