## PauliString that is a Pauli Operator
struct PauliString{Optype<:Integer}
    nqubits::Int
    operator::Optype
    coeff::Float64
end


function PauliString(nq, symbol::Symbol, qind::Int, coeff=1.0)
    inttype = getinttype(nq)
    temp_op = inttype(0)
    temp_op = setelement!(temp_op, qind, symboltoint(symbol))

    return PauliString(nq, symboltoint(temp_op), coeff)
end

function PauliString(nq, symbols::Vector{Symbol}, qinds::Vector{Int}, coeff=1.0)

    inttype = getinttype(nq)
    temp_op = inttype(0)
    for (op, qind) in zip(symbols, qinds)
        temp_op = setelement!(temp_op, qind, symboltoint(op))
    end

    return PauliString(nq, temp_op, coeff)
end

import Base.show
function show(io::IO, pstr::PauliString)
    pauli_string = inttostring(pstr.operator, pstr.nqubits)
    if length(pauli_string) > 20
        pauli_string = pauli_string[1:20] * "..."
    end
    print(io, "PauliString(nqubits: $(pstr.nqubits), $(round(pstr.coeff, sigdigits=5)) * $(pauli_string))")
end

## PauliSum that is a sum of PauliString
struct PauliSum{OpType<:Integer}
    nqubits::Int
    op_dict::Dict{OpType,Float64}
end

PauliSum(nq::Int) = PauliSum(nq, Dict{getinttype(nq),Float64}())

# PauliSum(nq::Int, coeff_type::DataType) = PauliSum(nq, Dict{getinttype(nq),coeff_type}())  # TODO: figure out how to use the numeric wrapper type :) 

function PauliSum(nq::Int, pstr::PauliString)
    # there is a possible downfall here nq != pauli_string.nqubits
    # this will convert to nq and consequently potentially truncate bits
    return PauliSum(nq, Dict{getinttype(nq),Float64}(pstr.operator => pstr.coeff))
end

function PauliSum(nq::Int, pstr::Vector{PauliString{T}}) where {T}
    # there is a possible downfall here nq != pauli_string.nqubits
    # this will convert to nq and consequently potentially truncate bits
    op_dict = Dict{getinttype(nq),Float64}(
        pauli_string.operator => pauli_string.coeff for pauli_string in pstr
    )
    return PauliSum(nq, op_dict)
end

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

## Adding to PauliSum
function add!(psum::PauliSum, pstr::PauliString)
    # there is a possible downfall here nq != pauli_string.nqubits
    # this will convert to nq and consequently potentially truncate bits

    if haskey(psum.op_dict, pstr.operator)
        psum.op_dict[pstr.operator] += pstr.coeff
    else
        psum.op_dict[pstr.operator] = pstr.coeff
    end

    return psum
end

function add!(psum::PauliSum, pstr::Vector{PauliString{T}}) where {T}
    # there is a possible downfall here nq != pauli_string.nqubits
    # this will convert to nq and consequently potentially truncate bits

    for pauli_string in pstr
        if haskey(psum.op_dict, pauli_string.operator)
            psum.op_dict[pauli_string.operator] += pauli_string.coeff
        else
            psum.op_dict[pauli_string.operator] = pauli_string.coeff
        end
    end

    return psum
end

function add!(psum, symbol, qind, coeff=1.0)
    return add!(psum, PauliString(psum.nqubits, symbol, qind, coeff))
end

function add!(psum, symbols::Vector{Symbol}, qinds::Vector{Int}, coeff=1.0)
    return add!(psum, PauliString(psum.nqubits, symbols, qinds, coeff))
end

## PathProperties # TODO: make compatible with PauliString and PauliSum
abstract type PathProperties end

mutable struct NumericPathProperties <: PathProperties
    coeff::Float64
    nsins::Int
    ncos::Int
    freq::Int
end

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