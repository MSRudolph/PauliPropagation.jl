###
##
# A file to conjugate/transpose our gates for Schrödinger picture propagation.
##
###

"""
    toschrodinger(circuit, params)

Convert a circuit and its parameters to the Schrödinger picture.
Calls `toschrodinger(gate [, param])` for each gate in the circuit
with an optional parameter if the gate is parametrized.
This function needs to be overloaded for custom gates.
Returns the converted circuit and parameters.
"""
function toschrodinger(circuit, params)
    params_iterator = Iterators.Stateful(params)
    schrodinger_circuit = eltype(circuit)[]
    schrodinger_params = eltype(params)[]
    for gate in circuit
        if gate isa ParametrizedGate
            param = popfirst!(params_iterator)
            gate_schrodinger, param_schrodinger = toschrodinger(gate, param)
            push!(schrodinger_circuit, gate_schrodinger)
            push!(schrodinger_params, param_schrodinger)
        else
            gate_schrodinger = toschrodinger(gate)
            push!(schrodinger_circuit, gate_schrodinger)
        end
    end
    return schrodinger_circuit, schrodinger_params
end

function toschrodinger(gate::G, args...) where G
    throw(error("Unkown how to define gate of type $G in the Schrodinger picture. 
    Please implement `toschrodinger(gate::G [, param])` for this gate type."))
end

"""
    toschrodinger(gate::PauliRotation, θ)

Method to transpose a `PauliRotation` gate for Schrödinger picture propagation.
This inverts the sign of `θ`, which can be seen by the conjugation of exp(-i*θ*P/2).
"""
function toschrodinger(gate::PauliRotation, θ)
    return gate, -θ
end

"""
    toschrodinger(gate::ImaginaryPauliRotation, τ)

Method to transpose a `ImaginaryPauliRotation` gate for Schrödinger picture propagation.
This does not change the parameter because exp(-τ*P/2) is the same under conjugation.
"""
function toschrodinger(gate::ImaginaryPauliRotation, τ)
    return gate, τ
end

"""
    toschrodinger(gate::CliffordGate)

Method to transpose a `CliffordGate` for Schrödinger picture propagation.
If not already registered, the transposed Clifford map is created via `transposecliffordmap()`
and stored in the global `clifford_map`. 
This Clifford gate is called `:(old_symbol)_transpose`, where `old_symbol` is the symbol of the original Clifford gate.
"""
function toschrodinger(gate::CliffordGate)
    transposed_symbol = Symbol(gate.symbol, :_transpose)

    if haskey(clifford_map, transposed_symbol)
        return gate
    end

    # register the transpose 
    lookup_map = clifford_map[gate.symbol]
    transposed_map = transposecliffordmap(lookup_map)
    clifford_map[transposed_symbol] = transposed_map

    return CliffordGate(transposed_symbol, gate.qinds)
end

"""
    toschrodinger(gate::PauliNoise, p)

Method to transpose a `PauliNoise` gate for Schrödinger picture propagation.
PauliNoise gates are symmetric under transposition and are thus unaffected.
"""
function toschrodinger(gate::PauliNoise, p)
    return gate, p
end


function toschrodinger(frozen_gate::FrozenGate)
    gate, param = frozen_gate.gate, frozen_gate.parameter
    gate_schrodinger, param_schrodinger = toschrodinger(gate, param)
    return freeze(gate_schrodinger, param_schrodinger)
end