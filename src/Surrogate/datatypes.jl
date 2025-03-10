"""
Abstract node type for the Pauli propagation Surrogate 
"""
abstract type CircuitNode end

"""
    EvalEndNode(pstr::Integer, coefficient::Real)

Node type for the Pauli strings in the observable to be backpropagated.
"""
@kwdef mutable struct EvalEndNode <: CircuitNode
    pstr::Int
    coefficient::Float64
    cummulative_value::Float64 = 0.0
    is_evaluated::Bool = false
end

"""
    EvalEndNode(pstr::Integer)

Constructor for `EvalEndNode` with a default coefficient of 1.0.
"""
EvalEndNode(pstr::Integer) = EvalEndNode(pstr, 1.0)

"""
    PauliRotationNode(parents::Vector{Union{EvalEndNode,PauliRotationNode}}, trig_inds::Vector{Int}, signs::Vector{Int}, param_idx::Int)

Surrogate graph node for a Pauli rotation gate.
"""
@kwdef mutable struct PauliRotationNode <: CircuitNode
    parents::Vector{Union{EvalEndNode,PauliRotationNode}}
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
Base.show(io::IO, node::EvalEndNode) = print(io, "$(typeof(node))(Pauli string=$(node.pstr), coefficient=$(node.coefficient))")

## PathProperties Type
"""
    NodePathProperties(node::CircuitNode, nsins::Int, ncos::Int, freq::Int)

Surrogate `PathProperties` type. Carries `CircuitNode`s instead of numerical coefficients.
"""
struct NodePathProperties <: PathProperties
    node::Union{EvalEndNode,PauliRotationNode}
    nsins::Int
    ncos::Int
    freq::Int
end

"""
Pretty print for PauliFreqTracker
"""
Base.show(io::IO, pth::NodePathProperties) = print(io, "NodePathProperties($(typeof(pth.node)), nsins=$(pth.nsins), ncos=$(pth.ncos), freq=$(pth.freq))")


"""
    NodePathProperties(node::CircuitNode)

One-argument constructor for `NodePathProperties`. Initializes `nsins`, `ncos`, and `freq` to 0.
"""
NodePathProperties(node::CircuitNode) = NodePathProperties(node, 0, 0, 0)


"""
    tonumber(path::NodePathProperties)

Get the cummulative coefficient of a `NodePathProperties` node.
This assumes that the surrogate has already been evaluated.
"""
tonumber(path::NodePathProperties) = path.node.cummulative_value


"""
    wrapcoefficients(pstr::PauliString, ::Type{NodePathProperties})

Wrap the coefficient of a `PauliString` into `NodePathProperties`.
"""
function wrapcoefficients(pstr::PauliString, ::Type{NodePathProperties})
    node = NodePathProperties(EvalEndNode(pstr.term, pstr.coeff, 0.0, false))
    return PauliString(pstr.nqubits, pstr.term, node)
end

"""
    wrapcoefficients(psum::PauliSum, ::Type{NodePathProperties})

Wrap the coefficients of a `PauliSum` into `NodePathProperties`.
"""
function wrapcoefficients(psum::PauliSum, ::Type{NodePathProperties})
    return PauliSum(psum.nqubits, Dict(pstr => NodePathProperties(EvalEndNode(pstr, coeff, 0.0, false)) for (pstr, coeff) in psum.terms))
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

"""
    set!(psum::Dict{TermType, NodePathProperties}, pstr::TermType, path::NodePathProperties)

In-place setting the coefficient of a Pauli string in a Pauli sum dictionary.
The type of the Pauli string needs to be the keytype=`TermType` of the dictionary, and the coefficient `coeff` needs to be the valuetype=`NodePathProperties`.
If the coefficient is 0, the Pauli string is removed from the dictionary.
"""
function set!(psum::Dict{TT,NodePathProperties}, pstr::TT, path::NodePathProperties) where {TT}
    psum[pstr] = path
    return psum
end