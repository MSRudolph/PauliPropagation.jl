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
            push!(circuit, PauliGate(:X, ii))

            # RZ
            push!(circuit, PauliGate(:Z, ii))

            # RX
            push!(circuit, PauliGate(:X, ii))
        end

        for pair in topology
            # CNOT or YY gate here
            push!(circuit, PauliGate([:Y, :Y], pair))
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
            push!(circuit, PauliGate(:Y, ii))

            # RZ
            push!(circuit, PauliGate(:Z, ii))

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

    zzlayer(circuit) = append!(circuit, (PauliGate([:Z, :Z], pair) for pair in topology))
    xlayer(circuit) = append!(circuit, (PauliGate(:X, ii) for ii in 1:nqubits))

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
            push!(circuit, PauliGate([:X, :X], pair))

            # YY
            push!(circuit, PauliGate([:Y, :Y], pair))

            # ZZ
            push!(circuit, PauliGate([:Z, :Z], pair))
        end
    end

    tofastgates!(circuit)
    return circuit
end

import Random: shuffle
"""
    trottercircuitandparams(hamil::PauliSum, order::Integer, nreps::Integer; random=true)

Create a circuit, and an accompanying list of gate parameters, which trotterises the unitary evolution operator (at time=1) of the given Hamiltonian.
These are given by the higher-order (as per `order`) symmetrized Suzuki-Trotter decomposition, using Childs-randomisation of each repetition (as per `nreps`), unless disabled by `random=false.`
The specified `order` must be 1, or an even integer.
This function returns tuple (circuit, params), where every gate in vector `circuit` is a FastPauliGate with corresponding angle given in `params`.
Because time appears as a simple prefactor of every PauliGate angle, simulating non-unity time involves merely scaling all `params` by the desired time.
For example, `propagate(circuit, pstr, 0.5 * params)` would simulate real-time evolution to time=`0.5`. 
"""
function trottercircuitandparams(hamil::PauliSum, order::Integer, nreps::Integer; random=true)

    # move each term into a separate PauliSum to indicate lack of commutation
    groups = [PauliSum(hamil.nqubits, Dict(str=>coeff)) for (str,coeff) in hamil]
    return trottercircuitandparams(groups, order, nreps; random=random)
end

"""
    trottercircuitandparams(commutinggroups::Vector{PauliSum}, order::Integer, nreps::Integer)

Creates a Trotter circuit in an identical fashion to `trottercircuitandparams(PauliSum, ...)`, but where the terms within each `PauliSum` are mapped to contiguous gates.
This enables specifying commuting groups of Pauli strings in a Hamiltonian, reducing the Trotter error.
Each `PauliSum` within `commutinggroups` is assumed to contain Pauli strings which all commute with one another.
As such, Trotterisation will not interweave the corresponding rotations of terms in distinct groups.  
"""
function trottercircuitandparams(commutinggroups::Vector{PauliSum{A,B}}, order::Integer, nreps::Integer; random=true) where {A,B}

    # validate inputs, just because it is easy to pass erroneous odd-order
    if (order < 1 || (order != 1 && order % 2 != 0))
        throw(ArgumentError("Argument 'order' must be a positive even integer, or one."))
    elseif (nreps < 1)
        throw(ArgumentError("Argument 'nreps' must be a positive integer."))
    end

    # we will build two lists; gates and their corresponding params...
    circuit::Vector{Gate} = []
    params::Vector{Number} = []

    # which are informed by orderings of the non-commuting PauliSums
    groups = [collect(group) for group in commutinggroups]
    
    # TODO:
    # below is not type-stable because this is merely circuit-preparation
    # code and not invoked in hot simulation loops. However, it is foreseeable
    # that users might repeatedly call this function to re-randomise their
    # Trotter circuits, warranting a performance improvement. So, ye who
    # understand typing - pls fix dis c:

    # inner-function to Suzuki-symmetrize a Hamiltonian into (str,coeff) pairs,
    # as per Hatano et al arXiv:math-ph/0506007
    function symmetrize(groups, phase, order)
        if (order == 1)
            return [[(str,coeff*phase) for (str,coeff) in group] for group in groups]
        elseif (order == 2)
            groups = symmetrize(groups, phase/2, 1)
            return [groups; reverse(groups)]
        end
        factor = 1/(4 - 4^(1/(order-1)))
        outer = symmetrize(groups, factor * phase, order-2)
        inner = symmetrize(groups, (1 - 4*factor) * phase, order-2)
        return [outer; outer; inner; outer; outer]
    end

    # produce trotter circuit by repeatedly symmetrizing a random hamil order,
    # as per Childs et al, Quantum 3, 182 (2019)
    for _ in 1:nreps
        groups = random ? shuffle(groups) : groups
        reordering = symmetrize(groups, 1/nreps, order)
        for group in reordering
            for (str, coeff) in group
                paulis, inds = getpaulisandinds(str)
                gate = PauliGate(paulis, inds) ################################## MUST UPDATE
                push!(circuit, gate)
                push!(params, coeff)
            end
        end
    end

    # cast all PauliGate to FastPauliGate
    tofastgates!(circuit)
    return (circuit, params)
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
    push!(circuit, PauliGate(:Z, pair[1]))
    push!(circuit, PauliGate(:X, pair[1]))
    push!(circuit, PauliGate(:Z, pair[1]))

    # arbitrary on q2
    push!(circuit, PauliGate(:Z, pair[2]))
    push!(circuit, PauliGate(:X, pair[2]))
    push!(circuit, PauliGate(:Z, pair[2]))

    # entanglers
    push!(circuit, PauliGate([:X, :X], pair))
    push!(circuit, PauliGate([:Y, :Y], pair))
    push!(circuit, PauliGate([:Z, :Z], pair))

    # arbitrary on q1
    push!(circuit, PauliGate(:Z, pair[1]))
    push!(circuit, PauliGate(:X, pair[1]))
    push!(circuit, PauliGate(:Z, pair[1]))

    # arbitrary on q2
    push!(circuit, PauliGate(:Z, pair[2]))
    push!(circuit, PauliGate(:X, pair[2]))
    push!(circuit, PauliGate(:Z, pair[2]))
end