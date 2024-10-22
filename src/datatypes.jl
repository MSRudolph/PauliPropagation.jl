## PauliString that is a Pauli Operator
# TODO: change OpType to Unsigned
struct PauliString{OpType<:Integer,CoeffType}
    nqubits::Int
    operator::OpType
    coeff::CoeffType
end

function PauliString(nq, symbol::Symbol, qind::Int, coeff=1.0)
    temp_op = symboltoint(nq, symbol, qind)
    coeff = convertcoefficients(coeff)
    return PauliString(nq, temp_op, coeff)
end

function PauliString(nq, symbols::Vector{Symbol}, qinds::Vector{Int}, coeff=1.0)
    temp_op = symboltoint(nq, symbols, qinds)
    coeff = convertcoefficients(coeff)
    return PauliString(nq, temp_op, coeff)
end

import Base.show
function show(io::IO, pstr::PauliString)
    pauli_string = inttostring(pstr.operator, pstr.nqubits)
    if length(pauli_string) > 20
        pauli_string = pauli_string[1:20] * "..."
    end
    if isa(pstr.coeff, Number)
        coeff_str = round(pstr.coeff, sigdigits=5)
    elseif isa(pstr.coeff, PathProperties)
        if isa(pstr.coeff, Number)
            coeff_str = "PathProperty($(round(pstr.coeff, sigdigits=5)))"
        else
            coeff_str = "PathProperty($(typeof(pstr.coeff)))"
        end
    else
        coeff_str = "($(typeof(pstr.coeff)))"
    end
    print(io, "PauliString(nqubits: $(pstr.nqubits), $(coeff_str) * $(pauli_string))")
end

## PauliSum that is a sum of PauliString
struct PauliSum{OpType<:Integer,CoeffType}
    nqubits::Int
    op_dict::Dict{OpType,CoeffType}
end

PauliSum(nq::Int) = PauliSum(nq, Dict{getinttype(nq),Float64}())

function PauliSum(nq::Int, sym_dict::Dict{Vector{Symbol},CoeffType}) where {CoeffType}
    """
    Construct a PauliSum from a dictionary of {symbols, coefficients}

    Args:
        nq: number of qubits
        sym_dict: dictionary of {symbols, coefficients}

    Returns:
        PauliSum
    """
    checknumberofqubits.(nq, keys(sym_dict))

    int_dict = Dict(symboltoint(k) => convertcoefficients(v) for (k, v) in sym_dict)

    return PauliSum(nq, int_dict)
end

PauliSum(pstr::PauliString) = PauliSum(pstr.nqubits, pstr)

PauliSum(psum::PauliSum) = psum

function PauliSum(nq::Int, pstr::PauliString{OpType,CoeffType}) where {OpType,CoeffType}
    checknumberofqubits(nq, pstr)

    return PauliSum(nq, Dict{OpType,CoeffType}(pstr.operator => pstr.coeff))
end

# get coefficients or terms in the PauliSum
function getcoeff(psum::PauliSum{OpType,CoeffType}, operator::Integer) where {OpType,CoeffType}
    return get(psum.op_dict, operator, CoeffType(0))
end

function getcoeff(psum::PauliSum{OpType,CoeffType}, operator::Symbol, qind::Int) where {OpType,CoeffType}
    return getcoeff(psum, symboltoint(psum.nqubits, operator, qind))
end

function getcoeff(psum::PauliSum{OpType,CoeffType}, operator::Vector{Symbol}, qinds) where {OpType,CoeffType}
    return getcoeff(psum, symboltoint(psum.nqubits, operator, qinds))
end

function getcoeff(psum::PauliSum{OpType,CoeffType}, operator::Vector{Symbol}) where {OpType,CoeffType}
    return getcoeff(psum, symboltoint(operator))
end

getpaulistrings(psum::PauliSum) = [PauliString(psum.nqubits, op, coeff) for (op, coeff) in psum.op_dict]

Base.copy(psum::PauliSum) = PauliSum(psum.nqubits, copy(psum.op_dict))

Base.iterate(psum::PauliSum, state=1) = iterate(psum.op_dict, state)


import Base.show
function show(io::IO, psum::PauliSum)
    if length(psum.op_dict) == 0
        dict_string = "(no operators)"
    else
        dict_string = getprettystr(psum.op_dict, psum.nqubits)
    end
    print(io, "PauliSum(nqubits: $(psum.nqubits), $dict_string)")
end

import Base.length
length(psum::PauliSum) = length(psum.op_dict)

import Base: ==
# Define equality for PauliSum
function ==(ps1::PauliSum, ps2::PauliSum)
    if ps1.nqubits != ps2.nqubits
        return false
    end

    return ps1.op_dict == ps2.op_dict
end

## Adding to PauliSum
function add(pobj1::Union{PauliSum,PauliString}, pobj2::Union{PauliSum,PauliString})
    # adds pobj2 to pobj1 in place
    checknumberofqubits(pobj1, pobj2)

    pobj1 = copy(pobj1) # or deepcopy?

    add!(pobj1, pobj1)

    return pobj1
end

function add!(psum::PauliSum, pstr::PauliString)
    checknumberofqubits(psum, pstr)

    if haskey(psum.op_dict, pstr.operator)
        psum.op_dict[pstr.operator] += pstr.coeff
    else
        psum.op_dict[pstr.operator] = pstr.coeff
    end

    return psum
end

function add!(psum1::PauliSum, psum2::PauliSum)
    checknumberofqubits(psum1, psum2)

    for (operator, coeff) in psum2.op_dict
        if haskey(psum1.op_dict, operator)
            psum1.op_dict[operator] += coeff
        else
            psum1.op_dict[operator] = coeff
        end
    end

    return psum1
end

function add!(psum, symbol, qind, coeff=1.0)
    return add!(psum, PauliString(psum.nqubits, symbol, qind, coeff))
end

function add!(psum, symbols::Vector{Symbol}, qinds::Vector{Int}, coeff=1.0)
    return add!(psum, PauliString(psum.nqubits, symbols, qinds, coeff))
end

## Subtracting from PauliSum
function subtract(psum1::PauliSum, psum2::PauliSum; precision::Float64=1e-10)
    # subtracts psum2 from psum1
    checknumberofqubits(psum1, psum2)
    psum1 = copy(psum1) # or deepcopy?

    subtract!(psum1, psum2, precision=precision)

    return psum1
end

function subtract(psum::PauliSum, pstr::PauliString; precision::Float64=1e-10)
    # subtracts psum2 from psum1
    checknumberofqubits(psum, pstr)

    subtract(psum, PauliSum(pstr), precision=precision)

    return psum1
end

function subtract!(psum1::PauliSum, psum2::PauliSum; precision::Float64=1e-10)
    checknumberofqubits(psum1, psum2)
    for (operator, coeff) in psum2.op_dict
        if haskey(psum1.op_dict, operator)
            psum1.op_dict[operator] -= coeff

            # Remove the operator if the resulting coefficient is small
            if abs(psum1.op_dict[operator]) < precision
                delete!(psum1.op_dict, operator)
            end

        else
            psum1.op_dict[operator] = -coeff
        end
    end

    return psum1
end

## Helper functions
function convertcoefficients(coeff)
    if isa(coeff, Integer)
        return Float64(coeff)
    elseif isa(coeff, Complex)
        return Complex{Float64}(coeff)
    else
        return coeff
    end
end

function checknumberofqubits(nq::Int, op::Union{PauliString,PauliSum})
    if nq != op.nqubits
        throw(
            ArgumentError(
                "Number of qubits ($(nq)) must equal number of qubits ($(op.nqubits)) in $(typeof(op))"
            )
        )
    end
end

function checknumberofqubits(nq::Int, op::Vector{Symbol})
    if nq != length(op)
        throw(
            ArgumentError(
                "Number of qubits ($(op1.nqubits)) must equal number of qubits ($(length(op))) in $(typeof(op))"
            )
        )
    end
end

function checknumberofqubits(op1::Union{PauliString,PauliSum}, op2::Union{PauliString,PauliSum})
    if op1.nqubits != op2.nqubits
        throw(
            ArgumentError(
                "Number of qubits ($(op1.nqubits)) in $(typeof(op1)) must equal number of qubits ($(op2.nqubits)) in $(typeof(op2))"
            )
        )
    end
end

## This type can be used to wrap coefficients and record custom properties
abstract type PathProperties end

## Specific PathProperties
mutable struct NumericPathProperties <: PathProperties
    coeff::Float64
    nsins::Int
    ncos::Int
    freq::Int
end

## This 1-argument constructor needs to be defined for any PathProperties type
NumericPathProperties(coeff) = NumericPathProperties(coeff, 0, 0, 0)

import Base: *
function *(pth::PathProperties, val::Number)
    pth.coeff *= val
    return pth
end
import Base: copy
function copy(path_properties::PathProperties)
    return typeof(path_properties)(path_properties.coeff, path_properties.nsins, path_properties.ncos, path_properties.freq)
end

Base.show(io::IO, pth::PathProperties) = print(io, "PathProperties($(typeof(pth.coeff)), nsins=$(pth.nsins), ncos=$(pth.ncos), freq=$(pth.freq))")

Base.show(io::IO, pth::NumericPathProperties) = print(io, "NumericPathProperties($(pth.coeff), nsins=$(pth.nsins), ncos=$(pth.ncos), freq=$(pth.freq))")


## Wrapping PauliString and PauliSum in PathProperties
function wrapcoefficients(pstr::PauliString, PathPropertiesType::Type{PP}) where {PP<:PathProperties}
    # the one-argument constructor of your PathProperties type must be defined
    return PauliString(pstr.nqubits, pstr.operator, PathPropertiesType(pstr.coeff))
end