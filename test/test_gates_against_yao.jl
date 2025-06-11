using Test
using LinearAlgebra
using Random
using Yao: X, Y, Z, H, Rx, Rz, Ry, chain, put, control, zero_state, expect, apply, rot, mat, matblock, swap, SWAP, time_evolve, kron
using PauliPropagation

# Gate Translation Functions

function _Rzz(θ)
    phases = [exp(-im * θ/2), exp(im * θ/2), exp(im * θ/2), exp(-im * θ/2)]
    return put(2, 1:2 => matblock(Diagonal(phases)))
end

function _clifford_to_yao(g::CliffordGate)
    symbol = g.symbol
    if symbol == :X
        return X
    elseif symbol == :Y
        return Y
    elseif symbol == :Z
        return Z
    elseif symbol == :H
        return H
    elseif symbol == :S
        return chain(1, Rz(π/2))
    elseif symbol == :SX
        return chain(1, Rx(π/2))
    elseif symbol == :SY
        return chain(1, Ry(π/2))
    elseif symbol == :CNOT
        return control(2, 1, 2=>X)
    elseif symbol == :CZ
        return control(2, 1, 2=>Z)
    elseif symbol == :SWAP
        return SWAP
    elseif symbol == :ZZpihalf
        return _Rzz(π/2)
    else
        error("Unsupported CliffordGate symbol: $symbol")
    end
end

function _build_yao_observable(symbols::Vector{Symbol}, qubits::Vector{Int}, nqubits::Int)
    length(symbols) == length(qubits) || throw(ArgumentError("Symbols and qubits must have same length"))
    all(1 .≤ qubits .≤ nqubits) || throw(ArgumentError("Qubit indices out of range"))
    blocks = map(zip(symbols, qubits)) do (sym, q)
        op = sym == :X ? X : sym == :Y ? Y : sym == :Z ? Z : 
             throw(ArgumentError("Unsupported Pauli symbol: $sym. Use :X, :Y, or :Z"))
        put(nqubits, q => op)
    end
    return length(blocks) == 1 ? only(blocks) : chain(blocks...)
end

function _yao_heisenberg_circuit(nqubits::Int, nsteps::Int, Jx::Float64, Jy::Float64, Jz::Float64, dt::Float64)
    yao_circ = chain(nqubits)
    for _ in 1:nsteps
        for i in 1:nqubits-1
            push!(yao_circ, put(nqubits, (i, i+1) => time_evolve(kron(X, X), -Jx*dt/2)))
            push!(yao_circ, put(nqubits, (i, i+1) => time_evolve(kron(Y, Y), -Jy*dt/2)))
            push!(yao_circ, put(nqubits, (i, i+1) => time_evolve(kron(Z, Z), -Jz*dt/2)))
        end
    end
    return yao_circ
end

function _yao_tfi_circuit(nqubits::Int, nsteps::Int, J::Float64, h::Float64, dt::Float64)
    yao_circ = chain(nqubits)
    for _ in 1:nsteps
        for i in 1:nqubits-1
            push!(yao_circ, put(nqubits, (i, i+1) => time_evolve(kron(Z, Z), -J*dt/2)))
        end
        for i in 1:nqubits
            push!(yao_circ, put(nqubits, i => Rx(-h*dt)))
        end
    end
    return yao_circ
end

function _register_inverse!(symbol::Symbol)
    inv_symbol = Symbol(string(symbol), "_inv")
    if !haskey(clifford_map, inv_symbol)
        original_map = clifford_map[symbol]
        clifford_map[inv_symbol] = transposecliffordmap(original_map)
    end
    return inv_symbol
end

function _invert_gates(gates::Vector{<:Any}, θs::Vector{Float64})
    rev_list = Any[]
    rev_θs   = Float64[]
    θ_index = length(θs)
    for gate in reverse(gates)
        if gate isa PauliRotation
            push!(rev_list, gate)
            push!(rev_θs, -θs[θ_index])
            θ_index -= 1
        elseif gate isa CliffordGate
            inv_sym = _register_inverse!(gate.symbol)
            push!(rev_list, CliffordGate(inv_sym, gate.qinds))
        else
            throw(ArgumentError("Cannot invert object of type $(typeof(gate))"))
        end
    end
    return rev_list, rev_θs
end

function _insert_gate!(gate_type, nqubits, rng, custom_gates, yao_ops, θs)
    q -> rand(rng, 1:nqubits)
    if gate_type == :H || gate_type == :X || gate_type == :Y || gate_type == :Z
        q = rand(rng, 1:nqubits)
        push!(custom_gates, CliffordGate(gate_type, [q]))
        g = gate_type == :H ? H : gate_type == :X ? X : gate_type == :Y ? Y : Z
        push!(yao_ops, put(nqubits, q => g))
    elseif gate_type == :RX || gate_type == :RY || gate_type == :RZ
        q = rand(rng, 1:nqubits)
        angle = rand(rng) * π/2
        axis = Symbol(string(gate_type)[end])
        push!(custom_gates, PauliRotation(axis, q))
        push!(θs, angle)
        rot = gate_type == :RX ? Rx(angle) : gate_type == :RY ? Ry(angle) : Rz(angle)
        push!(yao_ops, put(nqubits, q => rot))
    elseif gate_type == :CNOT && nqubits ≥ 2
        c, t = rand(rng, 1:nqubits), rand(rng, 1:nqubits)
        while t == c; t = rand(rng, 1:nqubits); end
        push!(custom_gates, CliffordGate(:CNOT, [c, t]))
        push!(yao_ops, control(nqubits, c, t => X))
    elseif gate_type == :SWAP && nqubits ≥ 2
        q1, q2 = rand(rng, 1:nqubits), rand(rng, 1:nqubits)
        while q2 == q1; q2 = rand(rng, 1:nqubits); end
        push!(custom_gates, CliffordGate(:SWAP, [q1, q2]))
        push!(yao_ops, swap(nqubits, q1, q2))
    elseif gate_type == :PauliRotation
        k = rand(rng, 1:min(3, nqubits))
        qs = sort(unique(rand(rng, 1:nqubits, k)))
        paulis = rand(rng, [:X, :Y, :Z], length(qs))
        angle = rand(rng) * π/2
        push!(custom_gates, PauliRotation(paulis, qs))
        push!(θs, angle)
        blocks = [pa == :X ? Rx(angle) : pa == :Y ? Ry(angle) : Rz(angle) for pa in paulis]
        yao_op = put(nqubits, Tuple(qs) => foldl(kron, blocks))
        push!(yao_ops, yao_op)
    end
end

const all_clifford_gates = collect(keys(PauliPropagation._default_clifford_map))
const single_obs = [:X, :Y, :Z]
const two_obs = [(:X, :X), (:X, :Y), (:X, :Z),
                 (:Y, :X), (:Y, :Y), (:Y, :Z),
                 (:Z, :X), (:Z, :Y), (:Z, :Z)]


# Test Clifford Gates on All Observables

@testset "Clifford Gate Propagation" begin
    for gate_symbol in all_clifford_gates
        nqubits = gate_symbol in (:CNOT, :CZ, :SWAP, :ZZpihalf) ? 2 : 1
        qinds = 1:nqubits
        @testset "$gate_symbol on Single Qubit Observables" begin
            yao_gate = _clifford_to_yao(CliffordGate(gate_symbol, qinds))
            for obs_symbol in single_obs
                qc = [CliffordGate(gate_symbol, qinds)]
                obs = PauliSum(nqubits)
                add!(obs, [obs_symbol], [1], 1)
                propagated = propagate(qc, obs)
                exp_val = overlapwithzero(propagated)
                state = zero_state(nqubits)
                evolved = apply(state, yao_gate)
                yao_obs = _build_yao_observable([obs_symbol], [1], nqubits)
                ref_val = real(expect(yao_obs, evolved))
                @test isapprox(exp_val, ref_val; atol=1e-10)
            end
        end
        nqubits == 2 && @testset "$gate_symbol on Two Qubit Observables" begin
            yao_gate = _clifford_to_yao(CliffordGate(gate_symbol, qinds))
            for (o1, o2) in two_obs
                qc = [CliffordGate(gate_symbol, qinds)]
                obs = PauliSum(2)
                add!(obs, [o1, o2], [1, 2], 1)
                propagated = propagate(qc, obs)
                exp_val = overlapwithzero(propagated)
                state = zero_state(2)
                evolved = apply(state, yao_gate)
                yao_obs = _build_yao_observable([o1, o2], [1, 2], 2)
                ref_val = real(expect(yao_obs, evolved))
                @test isapprox(exp_val, ref_val; atol=1e-10)
            end
        end
    end
end

# Test PauliRotation Gates

@testset "PauliRotation Gates" begin
    @testset "Basic Properties" begin
        @test tomatrix(PauliRotation(:Z, 1), 0) ≈ Matrix(I, 2, 2)
        @test tomatrix(PauliRotation(:Z, 1), π/2) ≈ [exp(-im*π/4) 0; 0 exp(im*π/4)]
        @test tomatrix(PauliRotation(:Z, 1), π) ≈ -im * mat(Z)
        @test tomatrix(PauliRotation(:X, 1), π) ≈ -im * mat(X)
        @test tomatrix(PauliRotation(:Y, 1), π) ≈ -im * mat(Y)
    end
    @testset "Against Yao Rotations" begin
        for (axis, yao_rot) in [(:X, Rx), (:Y, Ry), (:Z, Rz)]
            θ = randn()
            yao_gate = put(1, 1 => yao_rot(θ))
            pr_gate = PauliRotation(axis, 1)
            @test mat(yao_gate) ≈ tomatrix(pr_gate, θ)
        end
    end
    @testset "Multi-Qubit Rotations" begin
        for (symbols, op) in [([:X, :X], kron(X, X)),
                              ([:Y, :Y], kron(Y, Y)),
                              ([:Z, :Z], kron(Z, Z))]
            θ = randn()
            pr = PauliRotation(symbols, [1, 2])
            yao = put(2, (1,2) => time_evolve(op, θ/2))
            @test tomatrix(pr, θ) ≈ mat(yao)
        end
    end
end

# Test Random Circuits with PauliRotations and Cliffords

@testset "Randomized PauliRotation & Clifford Tests" begin
    rng = MersenneTwister(1234)
    for trial in 1:10
        nqubits = rand(rng, 1:3)
        depth = rand(rng, 3:6)
        custom_gates = Any[]
        yao_ops = Any[]
        θs = Float64[]
        for _ in 1:depth
            gate_type = rand(rng, [:H, :X, :Y, :Z, :RX, :RY, :RZ, :CNOT, :SWAP, :PauliRotation])
            _insert_gate!(gate_type, nqubits, rng, custom_gates, yao_ops, θs)
        end
        if rand(rng) < 0.5
            k = rand(rng, 1:nqubits)
            obs_qubits = sort(unique(rand(rng, 1:nqubits, k)))
            obs_symbols = rand(rng, [:X, :Y, :Z], length(obs_qubits))
        else
            obs_symbols = [:Z]
            obs_qubits = [rand(rng, 1:nqubits)]
        end
        @testset "Trial $trial (n=$nqubits, depth=$depth)" begin
            obs = PauliSum(nqubits)
            add!(obs, obs_symbols, obs_qubits, 1.0)
            propagated = propagate(custom_gates, obs, θs)
            custom_val = overlapwithzero(propagated)
            zero_st = zero_state(nqubits)
            evolved_state = apply(zero_st, chain(yao_ops...))
            yao_obs = _build_yao_observable(obs_symbols, obs_qubits, nqubits)
            yao_val = real(expect(yao_obs, evolved_state))
            @test isapprox(custom_val, yao_val; atol=1e-10)
            rev_gates, rev_θs = _invert_gates(custom_gates, θs)
            roundtrip_obs = propagate(rev_gates, propagated, rev_θs)
            orig_val_zero = overlapwithzero(obs)
            roundtrip_val_zero = overlapwithzero(roundtrip_obs)
            @test isapprox(roundtrip_val_zero, orig_val_zero; atol=1e-10)
        end
    end
end

# Test Model Hamiltonians

@testset "Transverse Field Ising Model" begin
    rng = MersenneTwister(42)
    for nqubits in [2, 3, 4]
        J = 1.0
        h = 0.5
        dt = 0.1
        nsteps = 3
        circuit = tfitrottercircuit(nqubits, nsteps)
        θs = Float64[]
        for _ in 1:nsteps
            for _ in 1:nqubits-1
                push!(θs, -J * dt)
            end
            for _ in 1:nqubits
                push!(θs, -h * dt)
            end
        end
        yao_circ = _yao_tfi_circuit(nqubits, nsteps, J, h, dt)
        state = zero_state(nqubits) |> yao_circ
        @testset "n=$nqubits" begin
            for q in 1:nqubits, p in [:X, :Y, :Z]
                obs = PauliSum(nqubits)
                add!(obs, [p], [q], 1.0)
                our_val = overlapwithzero(propagate(circuit, obs, θs))
                yao_obs = _build_yao_observable([p], [q], nqubits)
                yao_val = real(expect(yao_obs, state))
                @test isapprox(our_val, yao_val; atol=1e-10)
            end
            if nqubits ≥ 2
                for q1 in 1:nqubits-1
                    obs = PauliSum(nqubits)
                    add!(obs, [:Z, :Z], [q1, q1+1], 1.0)
                    our_val = overlapwithzero(propagate(circuit, obs, θs))
                    yao_obs = _build_yao_observable([:Z, :Z], [q1, q1+1], nqubits)
                    yao_val = real(expect(yao_obs, state))
                    @test isapprox(our_val, yao_val; atol=1e-10)
                end
            end
        end
    end
end

@testset "Heisenberg Model" begin
    rng = MersenneTwister(42)
    for nqubits in [2, 3]
        Jx, Jy, Jz = 0.8, 0.9, 1.0
        dt = 0.05
        nsteps = 2
        circuit = heisenbergtrottercircuit(nqubits, nsteps)
        θs = Float64[]
        for _ in 1:nsteps
            for _ in 1:nqubits-1
                push!(θs, -Jx * dt)
                push!(θs, -Jy * dt)
                push!(θs, -Jz * dt)
            end
        end
        yao_circ = _yao_heisenberg_circuit(nqubits, nsteps, Jx, Jy, Jz, dt)
        state = zero_state(nqubits) |> yao_circ
        @testset "n=$nqubits" begin
            for q in 1:nqubits, p in [:X, :Y, :Z]
                obs = PauliSum(nqubits)
                add!(obs, [p], [q], 1.0)
                our_val = overlapwithzero(propagate(circuit, obs, θs))
                yao_obs = _build_yao_observable([p], [q], nqubits)
                yao_val = real(expect(yao_obs, state))
                @test isapprox(our_val, yao_val; atol=1e-2)
            end
            if nqubits ≥ 2
                for (p1, p2) in [(:X, :X), (:Y, :Y), (:Z, :Z)]
                    obs = PauliSum(nqubits)
                    add!(obs, [p1, p2], [1, 2], 1.0)
                    our_val = overlapwithzero(propagate(circuit, obs, θs))
                    yao_obs = _build_yao_observable([p1, p2], [1, 2], nqubits)
                    yao_val = real(expect(yao_obs, state))
                    @test isapprox(our_val, yao_val; atol=1e-2)
                end
            end
        end
    end
end

# Integration Tests

@testset "Integration Tests for Circuits" begin
    XX = kron(X, X)
    YY = kron(Y, Y)
    ZZ = kron(Z, Z)
    XYZ = kron(X, kron(Y, Z))
    circuits = [
        (
            nqubits = 2,
            custom_gates = [PauliRotation(:X, 1)],
            yao_circuit  = chain(put(2, 1 => Rx(π/4))),
            obs          = ([:Z], [1])
          ),
        (
            nqubits = 2,
            custom_gates = [CliffordGate(:CNOT, [1, 2])],
            yao_circuit  = chain(control(2, 1, 2 => X)),
            obs          = ([:Z], [2])
        ),
        (
            nqubits = 2,
            custom_gates = [PauliRotation(:Z, 1)],
            yao_circuit  = chain(put(2, 1 => Rz(π/4))),
            obs          = ([:X], [1])
        ),
        (
            nqubits = 2,
            custom_gates = [PauliRotation([:X, :X], [1, 2])],
            yao_circuit  = chain(put(2, (1,2) => time_evolve(XX, π/8))),
            obs          = ([:Z, :Z], [1, 2])
        ),
        (
            nqubits = 2,
            custom_gates = [PauliRotation([:Y, :Y], [1, 2])],
            yao_circuit  = chain(put(2, (1,2) => time_evolve(YY, π/8))),
            obs          = ([:X, :X], [1, 2])
        ),
        (
            nqubits = 2,
            custom_gates = [PauliRotation([:Z, :Z], [1, 2])],
            yao_circuit  = chain(put(2, (1,2) => time_evolve(ZZ, π/8))),
            obs          = ([:Y, :Y], [1, 2])
        ),
        (
            nqubits = 3,
            custom_gates = [PauliRotation([:X, :Y, :Z], [1, 2, 3])],
            yao_circuit  = chain(put(3, (1,2,3) => time_evolve(XYZ, π/8))),
            obs          = ([:Z, :Y, :X], [1, 2, 3])
        ),
    ]
    for circuit in circuits
        @testset "nqubits=$(circuit.nqubits), obs=$(circuit.obs)" begin
            θs = fill(π/4, count(g -> g isa PauliRotation, circuit.custom_gates))
            obs = PauliSum(circuit.nqubits)
            obs_symbols, obs_qubits = circuit.obs
            add!(obs, obs_symbols, obs_qubits, 1.0)
            propagated = propagate(circuit.custom_gates, obs, θs)
            custom_val = overlapwithzero(propagated)
            zero_st = zero_state(circuit.nqubits)
            evolved_state = apply(zero_st, circuit.yao_circuit)
            yao_obs = _build_yao_observable(obs_symbols, obs_qubits, circuit.nqubits)
            yao_val = real(expect(yao_obs, evolved_state))
            @test isapprox(custom_val, yao_val; atol=1e-10)
            rev_gates, rev_θs = _invert_gates(circuit.custom_gates, θs)
            roundtrip_obs = propagate(rev_gates, propagated, rev_θs)
            @test roundtrip_obs == obs
        end
    end
end
