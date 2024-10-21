function countparameters(circuit::AbstractVector)
    n_params = 0
    for gate in circuit
        n_params += isa(gate, ParametrizedGate)
    end
    return n_params
end


function bricklayertopology(nq::Int; periodic=false)
    return bricklayertopology(1:nq; periodic=periodic)
end

function bricklayertopology(qindices; periodic=false)
    nq = length(qindices)
    if nq == 1
        return []
    elseif nq == 2
        return [(qindices[1], qindices[2])]
    else
        topo = []
        for ii in 1:2:nq-1
            push!(topo, (qindices[ii], qindices[ii+1]))
        end
        if periodic && (qindices[end] % 2 == 1)  # odd layer
            push!(topo, (qindices[end], qindices[1]))
        end
        for ii in 2:2:nq-1
            push!(topo, (qindices[ii], qindices[ii+1]))
        end
        if periodic && (qindices[end] % 2 == 0)  # even layer
            push!(topo, (qindices[end], qindices[1]))
        end

        return topo
    end
end

function get2dtopology(nx, ny)
    topology = []

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

    return topology

end

function get2dstaircasetopology(nx, ny)
    next_inds = [1]
    temp_inds = []

    topology = []
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

function hardwareefficientcircuit(n_qubits, n_layers; topology=nothing)
    circuit::Vector{Gate} = []

    if isnothing(topology)
        topology = [(ii, ii + 1) for ii in 1:n_qubits-1]
    end

    for jj in 1:n_layers
        for ii in 1:n_qubits
            # RX
            push!(circuit, PauliGate([:X], [ii]))

            # RZ
            push!(circuit, PauliGate([:Z], [ii]))

            # RX
            push!(circuit, PauliGate([:X], [ii]))
        end

        for pair in topology
            # CNOT or YY gate here
            push!(circuit, PauliGate([:Y, :Y], collect(pair)))
        end
    end

    tofastgates!(circuit)
    return circuit
end

function efficientsu2circuit(n_qubits, n_layers; topology=nothing)
    circuit::Vector{Gate} = []

    if isnothing(topology)
        topology = [(ii, ii + 1) for ii in 1:n_qubits-1]
    end

    for jj in 1:n_layers
        for ii in 1:n_qubits
            # RY
            push!(circuit, PauliGate([:Y], [ii]))  # TODO: make fast gates

            # RZ
            push!(circuit, PauliGate([:Z], [ii]))

        end

        for pair in topology
            # CNOT or YY gate here
            push!(circuit, CliffordGate(:CNOT, collect(pair)))
        end
    end

    tofastgates!(circuit)
    return circuit
end


function tfitrottercircuit(n_qubits, n_layers; topology=nothing, start_with_ZZ=true)
    circuit::Vector{Gate} = []

    if isnothing(topology)
        topology = [(ii, ii + 1) for ii in 1:n_qubits-1]
    end

    zzlayer(circuit) = append!(circuit, (PauliGate([:Z, :Z], collect(pair)) for pair in topology))
    xlayer(circuit) = append!(circuit, (PauliGate([:X], [ii]) for ii in 1:n_qubits))

    if start_with_ZZ
        zzlayer(circuit)
    end

    for _ in 1:n_layers-1
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

function heisenbergtrottercircuit(n_qubits, n_layers; topology=nothing)
    circuit::Vector{Gate} = []

    if isnothing(topology)
        topology = [(ii, ii + 1) for ii in 1:n_qubits-1]
    end

    for jj in 1:n_layers

        for pair in topology
            # XX
            push!(circuit, PauliGate([:X, :X], collect(pair)))

            # YY
            push!(circuit, PauliGate([:Y, :Y], collect(pair)))

            # ZZ
            push!(circuit, PauliGate([:Z, :Z], collect(pair)))
        end
    end

    tofastgates!(circuit)
    return circuit
end


function su4ansatz(n_qubits, n_layers; topology=nothing)
    circuit::Vector{Gate} = []

    if isnothing(topology)
        topology = [(ii, ii + 1) for ii in 1:n_qubits-1]
    end

    for nl in 1:n_layers
        for pair in topology
            appendSU4!(circuit, pair)
        end
    end

    tofastgates!(circuit)
    return circuit
end

function qcnnansatz(n_qubits; periodic=false)
    circuit::Vector{Gate} = []

    qselection = 1:n_qubits
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



function appendSU4!(circuit, pair)
    # arbitrary on q1
    push!(circuit, PauliGate([:Z], [pair[1]]))
    push!(circuit, PauliGate([:X], [pair[1]]))
    push!(circuit, PauliGate([:Z], [pair[1]]))

    # arbitrary on q2
    push!(circuit, PauliGate([:Z], [pair[2]]))
    push!(circuit, PauliGate([:X], [pair[2]]))
    push!(circuit, PauliGate([:Z], [pair[2]]))

    # entanglers
    push!(circuit, PauliGate([:X, :X], collect(pair)))
    push!(circuit, PauliGate([:Y, :Y], collect(pair)))
    push!(circuit, PauliGate([:Z, :Z], collect(pair)))

    # arbitrary on q1
    push!(circuit, PauliGate([:Z], [pair[1]]))
    push!(circuit, PauliGate([:X], [pair[1]]))
    push!(circuit, PauliGate([:Z], [pair[1]]))

    # arbitrary on q2
    push!(circuit, PauliGate([:Z], [pair[2]]))
    push!(circuit, PauliGate([:X], [pair[2]]))
    push!(circuit, PauliGate([:Z], [pair[2]]))
end