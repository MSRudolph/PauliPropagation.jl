abstract type Gate end

abstract type ParametrizedGate <: Gate end

abstract type StaticGate <: Gate end

include("pauligates.jl")
include("cliffordgates.jl")
include("noisechannels.jl")
