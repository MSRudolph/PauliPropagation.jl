### frozengates.jl
##
# A file for frozen gates. Frozen gates are wrappers around parametrized gates tied to a fixed parameter.
# This way a parameter does not need to be passed to `propagate`, because it is already attached to the gate.
# This is not the intended way of using parametrized gates if the parameters change of the gate parameters are to be optimized. 
##
###

"""
    FrozenGate(gate::ParametrizedGate, parameter::Number)

A `StaticGate` that wraps a `ParametrizedGate` with a fixed parameter.
These are used to fix the parameter of `ParametrizedGate` at the time of circuit construction.
This can be convenient but might exclude this parameter from being, e.g., differentiated by external libraries.
"""
struct FrozenGate{GateType<:ParametrizedGate,T<:Number} <: StaticGate
    gate::GateType
    parameter::T
end


function Base.show(io::IO, frozen_gate::FrozenGate)
    print(io, "FrozenGate($(frozen_gate.gate), parameter = $(round(frozen_gate.parameter, sigdigits=3)))")
end

"""
    freeze(gate::ParametrizedGate, parameter::Number)

Returns a `FrozenGate` wrapping the `gate` with the fixed `parameter`.
"""
function freeze(gate::ParametrizedGate, parameter::Number)
    return FrozenGate(gate, parameter)
end

"""
    freeze(gates::Vector{Gate}, parameters::Vector{Number})

Returns a vector of `Gate`s where `ParametrizedGate`s are frozen with their `parameters`.
"""
function freeze(gates, parameters)
    @assert countparameters(gates) == length(parameters)

    frozen_gates = Vector{Gate}(undef, length(gates))
    param_idx = 1
    for (i, gate) in enumerate(gates)
        if isa(gate, ParametrizedGate)
            frozen_gates[i] = FrozenGate(gate, parameters[param_idx])
            param_idx += 1
        else
            frozen_gates[i] = gate
        end
    end
    return frozen_gates
end