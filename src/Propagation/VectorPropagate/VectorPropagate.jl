module VectorPropagate

using PauliPropagation
using AcceleratedKernels
const AK = AcceleratedKernels
using Bits

include("datatypes.jl")
include("propagate.jl")
include("stateoverlap.jl")

export PropagationCache

end