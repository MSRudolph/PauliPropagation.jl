### Circuits/circuits.jl
##
# Functions for building circuits for unitary evolution.
# Circuits are generally vectors of `Gate` objects.
##
###

## Circuits corresponding to Trotter time evolution circuits
"""
    tfitrottercircuit(nqubits::Integer, nlayers::Integer; topology=bricklayertopology(nqubits), start_with_ZZ=true)

Create a circuit that corresponds to a Trotterization of the transverse-field Ising Hamiltonian. 
```math
H = ∑_{(i,j) ∈ topology} J_{ij} Z_i Z_j + ∑_{i=1 ... nqubits} h_i X_i
```
A topology can be specified as a list of pairs of qubit indices. If no topology is specified, a bricklayer topology is used.
If `start_with_ZZ` is set to `true`, the circuit starts with a layer of ZZ gates, else with a layer of X gates. This is relevant depending on the initial state.
"""
function tfitrottercircuit(nqubits::Integer, nlayers::Integer; topology=bricklayertopology(nqubits), start_with_ZZ=true)
    circuit::Vector{Gate} = []

    # the function after this expects a circuit with at least one layer and will always append something
    if nlayers == 0
        return circuit
    end

    if start_with_ZZ
        rzzlayer!(circuit, topology)
    end

    for _ in 1:nlayers-1
        rxlayer!(circuit, nqubits)
        rzzlayer!(circuit, topology)
    end

    rxlayer!(circuit, nqubits)

    if !start_with_ZZ
        rzzlayer!(circuit, topology)
    end

    return circuit
end

"""
    tiltedtfitrottercircuit(nqubits::Integer, nlayers::Integer; topology=bricklayertopology(nqubits))

Create a circuit that corresponds to a Trotterization of the transverse-field Ising Hamiltonian. 
```math
H = ∑_{(i,j) ∈ topology} J_{ij} Z_i Z_j 
  + ∑_{i=1 ... nqubits} h_i Z_i +  ∑_{i=1 ... nqubits} b_i X_i
```
A topology can be specified as a list of pairs of qubit indices. If no topology is specified, a bricklayer topology is used.
"""
function tiltedtfitrottercircuit(nqubits::Integer, layers::Integer; topology=bricklayertopology(nqubits))
    # TODO(YT): the ordering of the Trotter circuit layer will be important
    # for shallow circuits to minimize the Trotter error.
    circuit::Vector{Gate} = []

    for _ in 1:layers
        rzzlayer!(circuit, topology)
        rzlayer!(circuit, nqubits)
        rxlayer!(circuit, nqubits)
    end

    return circuit
end

"""
    heisenbergtrottercircuit(nqubits::Integer, nlayers::Integer; topology=bricklayertopology(nqubits))

Create a circuit that corresponds to a Trotterization of the Heisenberg Hamiltonian.
A topology can be specified as a list of pairs of qubit indices. If no topology is specified, a bricklayer topology is used.
Note that the gates are applied as layers of XX-YY-ZZ gates, not as layers of XX on all, then YY on all, then ZZ on all. On the bricklayer topology, these are equivalent.
"""
function heisenbergtrottercircuit(nqubits::Integer, nlayers::Integer; topology=bricklayertopology(nqubits))
    circuit::Vector{Gate} = []

    for jj in 1:nlayers
        rxxlayer!(circuit, topology)
        ryylayer!(circuit, topology)
        rzzlayer!(circuit, topology)
    end

    return circuit
end

## Generic circuits
"""
    hardwareefficientcircuit(nqubits::Integer, nlayers::Integer; topology=bricklayertopology(nqubits))

Create a hardware-efficient circuit consisting of layers of single-qubit X-Z-X Pauli gates and YY entangling gates.
A topology can be specified as a list of pairs of qubit indices. 
If no topology is specified, a bricklayer topology is used.
"""
function hardwareefficientcircuit(nqubits::Integer, nlayers::Integer; topology=bricklayertopology(nqubits))
    circuit::Vector{Gate} = []

    for _ in 1:nlayers
        for ii in 1:nqubits
            # RX
            push!(circuit, PauliRotation(:X, ii))

            # RZ
            push!(circuit, PauliRotation(:Z, ii))

            # RX
            push!(circuit, PauliRotation(:X, ii))
        end

        for pair in topology
            # CNOT or YY gate here
            push!(circuit, PauliRotation([:Y, :Y], pair))
        end
    end

    return circuit
end

"""
    efficientsu2circuit(nqubits::Integer, nlayers::Integer; topology=bricklayertopology(nqubits))

Create a hardware-efficient circuit consisting of layers of single-qubit Y-Z Pauli gates and CNOT entangling gates.
A topology can be specified as a list of pairs of qubit indices. If no topology is specified, a bricklayer topology is used.
"""
function efficientsu2circuit(nqubits::Integer, nlayers::Integer; topology=bricklayertopology(nqubits))
    # TODO: Technically the middle layers have Y-Z-Y Pauli gates, and the last Z-Y.
    circuit::Vector{Gate} = []

    for jj in 1:nlayers
        for ii in 1:nqubits
            # RY
            push!(circuit, PauliRotation(:Y, ii))

            # RZ
            push!(circuit, PauliRotation(:Z, ii))

        end

        for pair in topology
            # CNOT or YY gate here
            push!(circuit, CliffordGate(:CNOT, pair))
        end
    end

    return circuit
end

"""
    su4circuit(nqubits::Integer, nlayers::Integer; topology=bricklayertopology(nqubits))

Create a circuit that consists of layers of SU(4) gates on a given topology. 
SU(4) gates are decomposed via the KAK Decomposition into single-qubit Z-X-Z gates on each qubit, followed by XX-YY-ZZ gates and again single-qubit Z-X-Z gates, for a total of 15 Pauli gates.
A topology can be specified as a list of pairs of qubit indices. If no topology is specified, a bricklayer topology is used.
"""
function su4circuit(nqubits::Integer, nlayers::Integer; topology=bricklayertopology(nqubits))
    circuit::Vector{Gate} = []

    for nl in 1:nlayers
        for pair in topology
            appendSU4!(circuit, pair)
        end
    end

    return circuit
end

"""
    qcnncircuit(nqubits::Integer; periodic=false)

Create a Quantum Convolutional Neural Network (QCNN) circuit on `nqubits` qubits.
The topology for the circuit is created by creating bricklayer topologies on half the qubits every layer. The final qubits are qubit 1 and ~`nqubits/2`, which should be measured.
"""
function qcnncircuit(nqubits::Integer; periodic=false)
    circuit::Vector{Gate} = []

    qselection = 1:nqubits
    topology = []
    while length(qselection) > 1
        # @show qselection
        append!(topology, bricklayertopology(qselection; periodic=periodic))
        qselection = qselection[1:2:end]
    end

    for pair in topology
        appendSU4!(circuit, pair)
    end

    return circuit
end


"""
    appendSU4!(circuit, pair)

Append a layer of SU(4) gates to a circuit on a given pair of qubits.
The SU(4) gate is decomposed via the KAK Decomposition into single-qubit Z-X-Z gates on each qubit, followed by XX-YY-ZZ gates and again single-qubit Z-X-Z gates, for a total of 15 Pauli gates.
"""
function appendSU4!(circuit, pair)
    # arbitrary on q1
    push!(circuit, PauliRotation(:Z, pair[1]))
    push!(circuit, PauliRotation(:X, pair[1]))
    push!(circuit, PauliRotation(:Z, pair[1]))

    # arbitrary on q2
    push!(circuit, PauliRotation(:Z, pair[2]))
    push!(circuit, PauliRotation(:X, pair[2]))
    push!(circuit, PauliRotation(:Z, pair[2]))

    # entanglers
    push!(circuit, PauliRotation([:X, :X], pair))
    push!(circuit, PauliRotation([:Y, :Y], pair))
    push!(circuit, PauliRotation([:Z, :Z], pair))

    # arbitrary on q1
    push!(circuit, PauliRotation(:Z, pair[1]))
    push!(circuit, PauliRotation(:X, pair[1]))
    push!(circuit, PauliRotation(:Z, pair[1]))

    # arbitrary on q2
    push!(circuit, PauliRotation(:Z, pair[2]))
    push!(circuit, PauliRotation(:X, pair[2]))
    push!(circuit, PauliRotation(:Z, pair[2]))
end


"""
    rxlayer!(circuit, nqubits)

Append a layer of single-qubit `PauliRotation(:X, qind)` gates to the circuit for each pair in the topology.
"""
rxlayer!(circuit, nqubits) = append!(circuit, (PauliRotation(:X, ii) for ii in 1:nqubits))

"""
    rylayer!(circuit, nqubits)

Append a layer of single-qubit `PauliRotation(:Y, qind)` gates to the circuit for each pair in the topology.
"""
rylayer!(circuit, nqubits) = append!(circuit, (PauliRotation(:Y, ii) for ii in 1:nqubits))

"""
    rzlayer!(circuit, nqubits)

Append a layer of single-qubit `PauliRotation(:Z, qind)` gates to the circuit for each pair in the topology.
"""
rzlayer!(circuit, nqubits) = append!(circuit, (PauliRotation(:Z, ii) for ii in 1:nqubits))

"""
    rxxlayer!(circuit, topology)

Append a layer of two-qubit `PauliRotation([:X, :X], pair)` gates to the circuit for each pair in the topology.
"""
rxxlayer!(circuit, topology) = append!(circuit, (PauliRotation([:X, :X], pair) for pair in topology))

"""
    ryylayer!(circuit, topology)

Append a layer of two-qubit `PauliRotation([:Y, :Y], pair)` gates to the circuit for each pair in the topology.
"""
ryylayer!(circuit, topology) = append!(circuit, (PauliRotation([:Y, :Y], pair) for pair in topology))

"""
    rzzlayer!(circuit, topology)

Append a layer of two-qubit `PauliRotation([:Z, :Z], pair)` gates to the circuit for each pair in the topology.
"""
rzzlayer!(circuit, topology) = append!(circuit, (PauliRotation([:Z, :Z], pair) for pair in topology))