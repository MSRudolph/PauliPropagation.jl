using Bits
using BitIntegers

const PauliStringType = Union{UInt8,UInt16,UInt32,UInt64,UInt128,UInt256,BigInt,Int} # to be maintained when we adapt getinttype()
const PauliType = PauliStringType

include("datatypes.jl")
include("bitoperations.jl")
include("paulioperations.jl")
include("utils.jl")