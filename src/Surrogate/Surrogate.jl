### NOTE ###
# This is an experimental submodule of PauliPropagation.jl
# The Pauli propagation Surrogate works in the limited case of circuits consisting of Clifford gates and Pauli gates.
# The handling can be awkward and is not final.
# View this as legacy code.
# This code will be removed from this repo in the future.
############ 

include("datatypes.jl")
include("propagate.jl")
include("evaluate.jl")
