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
    coefficient::Float64
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
Pretty print for `CircuitNode`
"""
Base.show(io::IO, node::CircuitNode) = print(io, "$(typeof(node))($(length(node.parents)) parent(s), param_idx=$(node.param_idx))")
"""
Pretty print for `EvalEndNode`
"""
Base.show(io::IO, node::EvalEndNode) = print(io, "$(typeof(node))(Operator=$(node.operator), coefficient=$(node.coefficient))")

## PathProperties Type
"""
    NodePathProperties(coeff::CircuitNode, nsins::Int, ncos::Int, freq::Int)

Surrogate `PathProperties` type. Carries `CircuitNode`s instead of numerical coefficients.
"""
struct NodePathProperties <: PathProperties
    coeff::Union{EvalEndNode,PauliGateNode}
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
    getnumcoeff(val::NodePathProperties)

Get the cummulative coefficient of a `NodePathProperties` node.
This assumes that the surrogate has already been evaluated.
"""
getnumcoeff(node::NodePathProperties) = node.coeff.cummulative_value


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