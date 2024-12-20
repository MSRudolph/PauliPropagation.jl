using Bits
using BitIntegers

"""
    PauliStringType

A union type for the integer types used to represent Pauli strings.
Pauli strings are objects like X ⊗ Z ⊗ I ⊗ Y, where each term is a Pauli acting on a qubit.
"""
const PauliStringType = Union{UInt8,UInt16,UInt32,UInt64,UInt128,UInt256,BigInt,Int} # to be maintained when we adapt getinttype()

"""
    PauliType

A union type for the integer types used to represent Paulis.
Paulis, also known as Pauli operators, are objects like I, X, Y, Z acting on a single qubit.
"""
const PauliType = PauliStringType

include("datatypes.jl")
include("bitoperations.jl")
include("paulioperations.jl")
include("utils.jl")