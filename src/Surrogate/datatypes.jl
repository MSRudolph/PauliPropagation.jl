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
Base.show(io::IO, node::EvalEndNode) = print(io, "$(typeof(node))(Pauli string=$(inttostring(node.pstr)), coefficient=$(node.coefficient))")

## PathProperties Type
"""
    NodePathProperties(coeff::CircuitNode, nsins::Int, ncos::Int, freq::Int)

Surrogate `PathProperties` type. Carries `CircuitNode`s instead of numerical coefficients.
"""
struct NodePathProperties <: PathProperties
    coeff::Union{EvalEndNode,PauliRotationNode}
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
    numcoefftype(node::NodePathProperties)

Get the type of the numerical coefficient of a `NodePathProperties` object.
Returns the type of the `cummulative_value` field of the stored `CircuitNode`.
"""
numcoefftype(node::NodePathProperties) = typeof(tonumber(node))

"""
    numcoefftype(::Type{NodePathProperties})

Currently returns `Float64` as the type of the `cummulative_value` in `CircuitNode`s.
"""
numcoefftype(::Type{NodePathProperties}) = Float64



"""
    tonumber(val::NodePathProperties)

Get the cummulative coefficient of a `NodePathProperties` node.
This assumes that the surrogate has already been evaluated.
"""
tonumber(node::NodePathProperties) = node.coeff.cummulative_value


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
    return PauliSum(psum.nqubits, Dict(op => NodePathProperties(EvalEndNode(op, coeff, 0.0, false)) for (op, coeff) in psum.terms))
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
    set!(psum::Dict{TermType, NodePathProperties}, pstr::TermType, coeff::NodePathProperties)

In-place setting the coefficient of a Pauli string in a Pauli sum dictionary.
The type of the Pauli string needs to be the keytype=`TermType` of the dictionary, and the coefficient `coeff` needs to be the valuetype=`NodePathProperties`.
If the coefficient is 0, the Pauli string is removed from the dictionary.
"""
function set!(psum::Dict{TT,NodePathProperties}, pstr::TT, coeff::NodePathProperties) where {TT}
    psum[pstr] = coeff
    return psum
end