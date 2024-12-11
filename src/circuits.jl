#=
Functions for defining circuits for unitary evolution.
=#

#TODO: refactor this to add building blocks for circuit layers.

"""
    countparameters(circuit)

Utility function to count the number of gates of type `ParametrizedGate` in a circuit.
"""
function countparameters(circuit)
    nparams = 0
    for gate in circuit
        nparams += isa(gate, ParametrizedGate)
    end
    return nparams
end

## Topologies

"""
    bricklayertopology(nqubits::Integer; periodic=false)

Create the topology of a so-called 1D bricklayer circuit on `nqubits` qubits. It consists of two sublayers connecting odd-even and eve-odd qubit indices, respectively.
If `periodic` is set to `true`, the last qubit is connected to the first qubit.
"""
function bricklayertopology(nqubits::Integer; periodic=false)
    return bricklayertopology(1:nqubits; periodic=periodic)
end

"""
    bricklayertopology(qindices; periodic=false)

Create the topology of a so-called 1D bricklayer circuit on a subset of qubits indicated by `qindices`.
If `periodic` is set to `true`, the last qubit is connected to the first qubit.
"""
function bricklayertopology(qindices; periodic=false)
    nqubits = length(qindices)

    topology = Tuple{Int,Int}[]
    if nqubits == 1
        return topology
    elseif nqubits == 2
        push!(topology, (qindices[1], qindices[2]))
        return topology
    else
        for ii in 1:2:nqubits-1
            push!(topology, (qindices[ii], qindices[ii+1]))
        end
        if periodic && (qindices[end] % 2 == 1)  # odd layer
            push!(topology, (qindices[end], qindices[1]))
        end
        for ii in 2:2:nqubits-1
            push!(topology, (qindices[ii], qindices[ii+1]))
        end
        if periodic && (qindices[end] % 2 == 0)  # even layer
            push!(topology, (qindices[end], qindices[1]))
        end

        return topology
    end
end

"""
    staircasetopology(nqubits::Integer; periodic=false)

Create a 1D staircase topology on `nqubits` qubits. The qubits are connected in a staircase pattern, where qubit `i` is connected to qubit `i+1`.
If `periodic` is set to `true`, the last qubit is connected to the first qubit.
"""
function staircasetopology(nqubits::Integer; periodic=false)
    topology = [(ii, ii + 1) for ii in 1:nqubits-1]
    if periodic
        push!(topology, (nqubits, 1))
    end
    return topology
end

"""
    get2dtopology(nx::Integer, ny::Integer; periodic=false)

Create a 2D topology on a grid of `nx` by `ny` qubits. The order is one that works and may need to be adapted for specific purposes.
If `periodic` is set to `true`, the grid is connected periodically in both directions.
"""
function get2dtopology(nx::Integer, ny::Integer; periodic=false)
    topology = Tuple{Int,Int}[]

    for jj in 1:ny
        for ii in 1:nx

            if jj <= ny - 1
                push!(topology, ((jj - 1) * nx + ii, jj * nx + ii))
            end

            if ii + 1 <= nx
                push!(topology, ((jj - 1) * nx + ii, (jj - 1) * nx + ii + 1))
            end
        end
    end

    if periodic
        nq = nx * ny
        for ii in 1:nx
            push!(topology, (ii, nq - nx + ii))
        end


        for ii in 0:ny-1
            push!(topology, (ii * nx + 1, ii * nx + nx))
        end

        topology = [pair for pair in unique(topology) if pair[1] != pair[2]]
    end

    return topology

end

"""
    get2dstaircasetopology(nx::Integer, ny::Integer)

Create a 2D staircase topology on a grid of `nx` by `ny` qubits.
Mind the order of the topology, which forms a staircase spanning the grid -> in the Schrödinger picture <-. 
"""
function get2dstaircasetopology(nx::Integer, ny::Integer)
    next_inds = [1]
    temp_inds = []

    topology = Tuple{Int,Int}[]
    while length(next_inds) > 0
        for ind in next_inds
            if ind % nx != 0
                next_ind = ind + 1
                push!(topology, (ind, next_ind))
                push!(temp_inds, next_ind)
            end
            if ceil(Int, ind / nx) < ny
                next_ind = ind + nx
                push!(topology, (ind, next_ind))
                push!(temp_inds, next_ind)
            end
        end
        next_inds = temp_inds
        temp_inds = []

    end
    return unique(topology)
end

## Create circuits
"""
    hardwareefficientcircuit(nqubits::Integer, nlayers::Integer; topology=nothing)

Create a hardware-efficient circuit consisting of layers of single-qubit X-Z-X Pauli gates and YY entangling gates.
A topology can be specified as a list of pairs of qubit indices. 
If no topology is specified, a bricklayer topology is used.
"""
function hardwareefficientcircuit(nqubits::Integer, nlayers::Integer; topology=nothing)
    circuit::Vector{Gate} = []

    if isnothing(topology)
        topology = bricklayertopology(nqubits)
    end

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

    tofastgates!(circuit)
    return circuit
end

"""
    efficientsu2circuit(nqubits::Integer, nlayers::Integer; topology=nothing)

Create a hardware-efficient circuit consisting of layers of single-qubit Y-Z Pauli gates and CNOT entangling gates.
A topology can be specified as a list of pairs of qubit indices. If no topology is specified, a bricklayer topology is used.
"""
function efficientsu2circuit(nqubits::Integer, nlayers::Integer; topology=nothing)
    # TODO: Technically the middle layers have Y-Z-Y Pauli gates, and the last Z-Y.
    circuit::Vector{Gate} = []

    if isnothing(topology)
        topology = bricklayertopology(nqubits)
    end

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

    tofastgates!(circuit)
    return circuit
end

"""
    tfitrottercircuit(nqubits::Integer, nlayers::Integer; topology=nothing, start_with_ZZ=true)

Create a circuit that corresponds to a Trotterization of the transverse-field Ising Hamiltonian. 
A topology can be specified as a list of pairs of qubit indices. If no topology is specified, a bricklayer topology is used.
If `start_with_ZZ` is set to `true`, the circuit starts with a layer of ZZ gates, else with a layer of X gates. This is relevant depending on the initial state.
"""
function tfitrottercircuit(nqubits::Integer, nlayers::Integer; topology=nothing, start_with_ZZ=true)
    circuit::Vector{Gate} = []

    if isnothing(topology)
        topology = bricklayertopology(nqubits)
    end

    zzlayer(circuit) = append!(circuit, (PauliRotation([:Z, :Z], pair) for pair in topology))
    xlayer(circuit) = append!(circuit, (PauliRotation(:X, ii) for ii in 1:nqubits))

    if start_with_ZZ
        zzlayer(circuit)
    end

    for _ in 1:nlayers-1
        xlayer(circuit)
        zzlayer(circuit)
    end

    xlayer(circuit)

    if !start_with_ZZ
        zzlayer(circuit)
    end

    tofastgates!(circuit)
    return circuit
end

"""
    tiltedtfitrottercircuit(n_qubits, n_layers; topology=nothing)

Returns a Trottterized circuit for the tilted transverse field Ising model.
H = Sum_{(i, i+1) in topology} Z_i Z_{i+1} 
  + Sum_{i=1}^{n_qubits} Z_i + Sum_{i=1}^{n_qubits} X_i

# Arguments
- `n_qubits::Integer`: The number of qubits in the circuit.
- `n_layers::Integer`: The number of Trotter steps to perform.
- `topology=nothing`: The topology of the qubits in the circuit. 
    Default (nothing): A linear chain.

# Returns
The Trottterized circuit as a vector of Gate.
"""
function tiltedtfitrottercircuit(n_qubits, n_layers; topology=nothing)
    # TODO(YT): the ordering of the Trotter circuit layer will be important
    # for shallow circuits to minimize the Trotter error.
    circuit::Vector{Gate} = []

    if isnothing(topology)
        topology = [(ii, ii + 1) for ii in 1:n_qubits-1]
    end

    zzlayer(circuit) = append!(
        circuit, (PauliRotation([:Z, :Z], pair) for pair in topology)
    )
    zlayer(circuit) = append!(
        circuit, (PauliRotation([:Z], [ii]) for ii in 1:n_qubits)
    )
    xlayer(circuit) = append!(
        circuit, (PauliRotation([:X], [ii]) for ii in 1:n_qubits)
    )

    for _ in 1:n_layers
        zzlayer(circuit)
        zlayer(circuit)
        xlayer(circuit)
    end

    tofastgates!(circuit)

    return circuit
end

"""
    heisenbergtrottercircuit(nqubits::Integer, nlayers::Integer; topology=nothing)

Create a circuit that corresponds to a Trotterization of the Heisenberg Hamiltonian.
A topology can be specified as a list of pairs of qubit indices. If no topology is specified, a bricklayer topology is used.
Note that the gates are applied as layers of XX-YY-ZZ gates, not as layers of XX on all, then YY on all, then ZZ on all. On the bricklayer topology, these are equivalent.
"""
function heisenbergtrottercircuit(nqubits::Integer, nlayers::Integer; topology=nothing)
    circuit::Vector{Gate} = []

    if isnothing(topology)
        topology = bricklayertopology(nqubits)
    end

    for jj in 1:nlayers

        for pair in topology
            # XX
            push!(circuit, PauliRotation([:X, :X], pair))

            # YY
            push!(circuit, PauliRotation([:Y, :Y], pair))

            # ZZ
            push!(circuit, PauliRotation([:Z, :Z], pair))
        end
    end

    tofastgates!(circuit)
    return circuit
end

"""
    su4ansatz(nqubits::Integer, nlayers::Integer; topology=nothing)

Create a circuit that consists of layers of SU(4) gates on a given topology. 
SU(4) gates are decomposed via the KAK Decomposition into single-qubit Z-X-Z gates on each qubit, followed by XX-YY-ZZ gates and again single-qubit Z-X-Z gates, for a total of 15 Pauli gates.
A topology can be specified as a list of pairs of qubit indices. If no topology is specified, a bricklayer topology is used.
"""
function su4ansatz(nqubits::Integer, nlayers::Integer; topology=nothing)
    circuit::Vector{Gate} = []

    if isnothing(topology)
        topology = bricklayertopology(nqubits)
    end

    for nl in 1:nlayers
        for pair in topology
            appendSU4!(circuit, pair)
        end
    end

    tofastgates!(circuit)
    return circuit
end

"""
    qcnnansatz(nqubits::Integer; periodic=false)

Create a Quantum Convolutional Neural Network (QCNN) ansatz on `nqubits` qubits.
The topology for the ansatz is created by creating bricklayer topologies on half the qubits every layer. The final qubits are qubit 1 and ~`nqubits/2`, which should be measured.
"""
function qcnnansatz(nqubits::Integer; periodic=false)
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

    tofastgates!(circuit)
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


const ibmeagletopology = [
    (1, 2), (1, 15), (2, 3), (3, 4), (4, 5), (5, 6), (5, 16), (6, 7), (7, 8), (8, 9), (9, 10), (9, 17),
    (10, 11), (11, 12), (12, 13), (13, 14), (13, 18), (15, 19), (16, 23), (17, 27), (18, 31), (19, 20),
    (20, 21), (21, 22), (21, 34), (22, 23), (23, 24), (24, 25), (25, 26), (25, 35), (26, 27), (27, 28),
    (28, 29), (29, 30), (29, 36), (30, 31), (31, 32), (32, 33), (33, 37), (34, 40), (35, 44), (36, 48),
    (37, 52), (38, 39), (38, 53), (39, 40), (40, 41), (41, 42), (42, 43), (42, 54), (43, 44), (44, 45),
    (45, 46), (46, 47), (46, 55), (47, 48), (48, 49), (49, 50), (50, 51), (50, 56), (51, 52), (53, 57),
    (54, 61), (55, 65), (56, 69), (57, 58), (58, 59), (59, 60), (59, 72), (60, 61), (61, 62), (62, 63),
    (63, 64), (63, 73), (64, 65), (65, 66), (66, 67), (67, 68), (67, 74), (68, 69), (69, 70), (70, 71),
    (71, 75), (72, 78), (73, 82), (74, 86), (75, 90), (76, 77), (76, 91), (77, 78), (78, 79), (79, 80),
    (80, 81), (80, 92), (81, 82), (82, 83), (83, 84), (84, 85), (84, 93), (85, 86), (86, 87), (87, 88),
    (88, 89), (88, 94), (89, 90), (91, 95), (92, 99), (93, 103), (94, 107), (95, 96), (96, 97), (97, 98),
    (97, 110), (98, 99), (99, 100), (100, 101), (101, 102), (101, 111), (102, 103), (103, 104), (104, 105),
    (105, 106), (105, 112), (106, 107), (107, 108), (108, 109), (109, 113), (110, 115), (111, 119), (112, 123),
    (113, 127), (114, 115), (115, 116), (116, 117), (117, 118), (118, 119), (119, 120), (120, 121), (121, 122),
    (122, 123), (123, 124), (124, 125), (125, 126), (126, 127)
]
