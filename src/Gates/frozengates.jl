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

import Base.show
function show(io::IO, frozen_gate::FrozenGate)
    print(io, "FrozenGate($(frozen_gate.gate), θ = $(round(frozen_gate.parameter, sigdigits=3)))")
end

"""
    freeze(gate::ParametrizedGate, parameter::Number)

Returns a `FrozenGate` wrapping the `gate` with the fixed `parameter`.
"""
function freeze(gate::ParametrizedGate, parameter::Number)
    return FrozenGate(gate, parameter)
end

"""
    freeze(gates, parameters)

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

"""
Get the maximum qubit index that the frozen gate acts on.
"""
_getmaxqubits(frozen_gate::FrozenGate) = _getmaxqubits(frozen_gate.gate)

"""
    apply(frozen_gate::FrozenGate, pstr, theta, args...)

Apply a `FrozenGate` to a Pauli string `pstr` with the parameter `FrozenGate.parameter`.
The passed `theta` will be ignored.
"""
function apply(frozen_gate::FrozenGate, pstr, theta, coefficient=1.0; kwargs...)
    return apply(frozen_gate.gate, pstr, frozen_gate.parameter, coefficient; kwargs...)
end

"""
Get the maxium number of qubits a gate acts on by calling this function on the wrapped gate.
"""
_maxqubits(gate::FrozenGate) = _maxqubits(gate.gate)

