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

"""
    apply(gate::StaticGate, pstr, theta, coefficient=1.0)

Calling apply on a `StaticGate` will dispatch to a 3-argument apply function without the paramter `theta`.
If a 4-argument apply function is defined for a concrete type, it will still dispatch to that one.
"""
apply(gate::SG, pstr, theta, coefficient=1.0; kwargs...) where {SG<:StaticGate} = apply(gate, pstr, coefficient; kwargs...)


include("frozengates.jl")
include("paulirotations.jl")
include("cliffordgates.jl")
include("noisechannels.jl")
include("miscgates.jl")
