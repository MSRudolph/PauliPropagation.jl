### PathProperties.jl
##
# PathProperties is a type that you can use to keep track of certain properties during propagation
# Currently provide only one concrete type PauliFreqTracker that keeps track of the number of sin and cos factors applied via a PauliRotation gate.
##
###

include("abstracttype.jl")
include("paulifreqtracker.jl")
