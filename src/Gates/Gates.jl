"""
Abstract type for gates. 
"""
abstract type Gate end

"""
Abstract type for parametrized gates.
"""
abstract type ParametrizedGate <: Gate end

"""
Abstract type for static gates are not parametrized.
"""
abstract type StaticGate <: Gate end

include("pauligates.jl")
include("cliffordgates.jl")
include("noisechannels.jl")
