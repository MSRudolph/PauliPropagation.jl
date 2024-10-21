## PauliString that is a Pauli Operator
# TODO: change OpType to Unsigned
struct PauliString{OpType<:Integer,CoeffType}
    nqubits::Int
    operator::OpType
    coeff::CoeffType
end

# TODO: Do things go wrong if people provide integer coefficients?
# Add a wrapper type `WrappedPauliString` that can hold `PathProperties` vs keep them separate as currently done?
function PauliString(nq, symbol::Symbol, qind::Int, coeff=1.0)
    inttype = getinttype(nq)
    temp_op = inttype(0)
    temp_op = setelement!(temp_op, qind, symboltoint(symbol))
    #TODO: find a better way to handle integer coefficients
    if isa(coeff, Int)  # Convert it to Float64 to have a stable type
        coeff = Float64(coeff) 
    end

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

function PauliSum(nq::Int, sym_dict::Dict{Vector{Symbol}, Float64})
    """
    Construct a PauliSum from a dictionary of {symbols, coefficients}

    Args:
        nq: number of qubits
        sym_dict: dictionary of {symbols, coefficients}

    Returns:
        PauliSum
    """
    int_dict = Dict(symboltoint(k) => v for (k, v) in sym_dict)

    return PauliSum(nq, int_dict)
end

function PauliSum(nq::Int, pstr::PauliString)
    if nq != pstr.nqubits
        throw(
            ArgumentError(
                "Number of qubits $nq must equal Pauli string $(pstr.nqubits)"
            )
        )
    end

    return PauliSum(nq, Dict{getinttype(nq),typeof(pstr.coeff)}(pstr.operator => pstr.coeff))
end

function PauliSum(nq::Int, pstr::Vector{PauliString{T1,T2}}) where {T1,T2}
    # there is a possible downfall here nq != pauli_string.nqubits
    # this will convert to nq and consequently potentially truncate bits
    op_dict = Dict{getinttype(nq),Float64}(
        pauli_string.operator => pauli_string.coeff for pauli_string in pstr
    )
    return PauliSum(nq, op_dict)
end

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

function add(psum1::PauliSum, psum2::PauliSum)
    # adds psum2 to psum1. Non-mutating version only defined between two PauliSum objects

    psum1 = copy(psum1) # or deepcopy?

    add!(psum1, psum2)

    return psum1
end

function add!(psum1::PauliSum, psum2::PauliSum)
    # there is a possible downfall here nq != pauli_string.nqubits
    # this will use nq from psum1

    for (operator, coeff) in psum2.op_dict
        if haskey(psum1.op_dict, operator)
            psum1.op_dict[operator] += coeff
        else
            psum1.op_dict[operator] = coeff
        end
    end

    return psum1
end

function add!(psum::PauliSum, pstr::Vector{PauliString{T1,T2}}) where {T1,T2}
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

#TODO: add support to subtract PauliString from PauliSum
function subtract(psum1::PauliSum, psum2::PauliSum; precision::Float64=1e-10)
    # subtracts psum2 from psum1

    psum1 = copy(psum1) # or deepcopy?

    subtract!(psum1, psum2, precision=precision)

    return psum1
end

function subtract!(psum1::PauliSum, psum2::PauliSum; precision::Float64=1e-10)
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