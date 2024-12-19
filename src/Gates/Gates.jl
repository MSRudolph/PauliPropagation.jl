### Gates.jl
##
# The top level file for gates.
# Gates are defined as structs that subtype either `ParametrizedGate` or `StaticGate`.
# How the gates act is defined in the `Propagation` module.
##
###

"""
Abstract type for gates. 
"""
abstract type Gate end

"""
Abstract type for parametrized gates.
"""
abstract type ParametrizedGate <: Gate end

"""
Abstract type for static gates that are not parametrized.
"""
abstract type StaticGate <: Gate end


include("frozengates.jl")
include("paulirotations.jl")
include("cliffordgates.jl")
include("noisechannels.jl")
include("miscgates.jl")
