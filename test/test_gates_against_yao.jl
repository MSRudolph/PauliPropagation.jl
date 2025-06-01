using Test
using LinearAlgebra
using Yao
using Random
using Yao: X, Y, Z, H, chain, put, control, zero_state, expect, apply

function Rzz(θ)
    phases = [exp(-im * θ/2), exp(im * θ/2), exp(im * θ/2), exp(-im * θ/2)]
    return put(2, 1:2 => matblock(Diagonal(phases)))
end

function clifford_to_yao(g::CliffordGate)
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
        return Rzz(π/2)
    else
        error("Unsupported CliffordGate symbol: $symbol")
    end
end

const all_gates = collect(keys(PauliPropagation._default_clifford_map))
const obs_symbols = [:X, :Y, :Z]

@testset "Observable Propagation Check via Yao.jl" begin
    for gate_symbol in all_gates
        nqubits = gate_symbol == :CNOT || gate_symbol == :CZ || gate_symbol == :SWAP || gate_symbol == :ZZpihalf ? 2 : 1
        for obs_symbol in obs_symbols
            @testset "Gate $gate_symbol propagating $obs_symbol" begin
                qc = [CliffordGate(gate_symbol, 1:nqubits)]
                obs = PauliSum(nqubits)
                add!(obs, [obs_symbol], [1], 1)
                propagated = propagate(qc, obs)
                exp_val_custom = overlapwithzero(propagated)
                state = zero_state(nqubits)
                yao_chain = chain(nqubits, [clifford_to_yao(g) for g in qc])
                evolved_state = Yao.apply(state, yao_chain)
                single_obs = obs_symbol == :X ? X : obs_symbol == :Y ? Y : Z
                yao_obs = put(nqubits, 1 => single_obs)
                ref_val = real(expect(yao_obs, evolved_state))
                @test isapprox(exp_val_custom, ref_val, atol=1e-10)
            end
        end
    end

    rot_gates = [
        (:X, Rx),
        (:Y, Ry),
        (:Z, Rz)
    ]
    for (axis, rot_gate_constructor) in rot_gates
        for test_i in 1:50
            θ = 2π * rand()
            pr_gate = PauliRotation(axis, 1)
            yao_gate = put(1 => rot_gate_constructor(θ))
            @test mat(chain(1, yao_gate)) ≈ tomatrix(pr_gate, θ) atol=1e-12
            for obs_symbol in obs_symbols
                @testset "Rotation Gate R$axis($θ) propagating $obs_symbol" begin
                    qc = [PauliRotation(axis, 1)]
                    obs = PauliSum(1)
                    add!(obs, [obs_symbol], [1], 1)
                    propagated = propagate(qc, obs, [θ])
                    exp_val_custom = overlapwithzero(propagated)
                    state = zero_state(1)
                    evolved_state = Yao.apply(state, yao_gate)
                    single_obs = obs_symbol == :X ? X : obs_symbol == :Y ? Y : Z
                    yao_obs = put(1, 1 => single_obs)
                    ref_val = real(expect(yao_obs, evolved_state))
                    @test isapprox(exp_val_custom, ref_val, atol=1e-10)
                end
            end
        end
    end
end

@testset "Combined Gates Propagation Test" begin
    nqubits = 2
    θ = π/4
    obs_symbols = [:X, :Y, :Z]
    qc = [
        CliffordGate(:H, [1]),
        PauliRotation(:Z, 1),
        CliffordGate(:CNOT, [1, 2])
    ]
    yao_gates = [
        put(nqubits, 1 => H),
        put(nqubits, 1 => Rz(θ)),
        control(nqubits, 1, 2 => X)
    ]
    yao_chain = chain(nqubits, yao_gates)
    state = zero_state(nqubits)
    for obs_symbol in obs_symbols
        @testset "Combined gates propagating $obs_symbol on qubit 1" begin
            obs = PauliSum(nqubits)
            add!(obs, [obs_symbol], [1], 1)
            propagated = propagate(qc, obs, [θ])
            exp_val_custom = overlapwithzero(propagated)
            evolved_state = Yao.apply(state, yao_chain)
            single_obs = obs_symbol == :X ? X : obs_symbol == :Y ? Y : Z
            yao_obs = put(nqubits, 1 => single_obs)
            ref_val = real(expect(yao_obs, evolved_state))
            @test isapprox(exp_val_custom, ref_val, atol=1e-10)
        end
    end
end

@testset "Pauli Rotation Gates" begin
    @testset "Single Qubit Rotations" begin
        @test tomatrix(PauliRotation(:Z, 1), 0) ≈ Matrix(I, 2, 2)
        @test tomatrix(PauliRotation(:Z, 1), π/2) ≈ [exp(-im*π/4) 0; 0 exp(im*π/4)]
        @test tomatrix(PauliRotation(:Z, 1), π) ≈ -im * mat(Z)
        @test tomatrix(PauliRotation(:X, 1), π) ≈ -im * mat(X)
        @test tomatrix(PauliRotation(:X, 1), π/2) ≈ [cos(π/4) -im*sin(π/4); -im*sin(π/4) cos(π/4)]
        @test tomatrix(PauliRotation(:Y, 1), π) ≈ -im * mat(Y)
        @test tomatrix(PauliRotation(:Y, 1), π/2) ≈ [cos(π/4) -sin(π/4); sin(π/4) cos(π/4)]
    end

    for gate in [(:X, Rx), (:Y, Ry), (:Z, Rz)]
        θ = randn()
        yao_gate = chain(1, put(1=>gate[2](θ)))
        pr_gate = PauliRotation(gate[1], 1)
        @test mat(yao_gate) ≈ tomatrix(pr_gate, θ)
    end
end

const SX = rot(X, π/2)
const SY = rot(Y, π/2)
const S  = Rz(π/2)

function invert_gates(gates::Vector{<:Any}, θs::Vector{Float64})
    rev_gates = reverse(gates)
    rev_θs = Float64[]
    θ_index = length(θs)
    for gate in rev_gates
        if gate isa PauliRotation
            push!(rev_θs, -θs[θ_index])
            θ_index -= 1
        end
    end

    return rev_gates, rev_θs
end

function build_yao_observable(symbols::Vector{Symbol}, qubits::Vector{Int}, nqubits::Int)
    blocks = []
    for (sym, q) in zip(symbols, qubits)
        op = sym == :X ? X : sym == :Y ? Y : sym == :Z ? Z : error("Unsupported symbol: $sym")
        push!(blocks, put(nqubits, q => op))
    end
    return length(blocks) == 1 ? blocks[1] : chain(blocks...)
end

@testset "Integration Tests for Full Circuits" begin
    circuits = [
        (
            name = "Simple H + CNOT + Z rotation",
            nqubits = 2,
            custom_gates = [
                CliffordGate(:H, [1]),
                CliffordGate(:CNOT, [1, 2]),
                PauliRotation(:Z, 2)
            ],
            yao_circuit = chain(
                put(2, 1 => H),
                control(2, 1, 2 => X),
                put(2, 2 => Rz(π/4))
            ),
            obs = ([:X], [1])
        ),
        (
            name = "Entangled ZZ circuit",
            nqubits = 2,
            custom_gates = [
                CliffordGate(:H, [1]),
                CliffordGate(:H, [2]),
                CliffordGate(:ZZpihalf, [1, 2])
            ],
            yao_circuit = chain(
                put(2, 1 => H),
                put(2, 2 => H),
                control(2, 1, 2 => Rz(π/2))
            ),
            obs = ([:Z], [1])
        ),
        (
            name = "SWAP + Y Rotation",
            nqubits = 2,
            custom_gates = [
                CliffordGate(:H, [1]),
                CliffordGate(:SWAP, [1, 2]),
                PauliRotation(:Y, 2)
            ],
            yao_circuit = chain(
                put(2, 1 => H),
                swap(2, 1, 2),
                put(2, 2 => Ry(π/4))
            ),
            obs = ([:Y], [1])
        ),
        (
            name = "3-qubit entangling circuit",
            nqubits = 3,
            custom_gates = [
                CliffordGate(:H, [1]),
                CliffordGate(:CNOT, [1, 2]),
                CliffordGate(:CNOT, [2, 3])
            ],
            yao_circuit = chain(
                put(3, 1 => H),
                control(3, 1, 2 => X),
                control(3, 2, 3 => X)
            ),
            obs = ([:X], [3])
        ),
        (
            name = "Mixed Clifford + Parametric",
            nqubits = 2,
            custom_gates = [
                CliffordGate(:H, [1]),
                PauliRotation(:X, 1),
                CliffordGate(:CNOT, [1, 2]),
                PauliRotation(:Y, 2)
            ],
            yao_circuit = chain(
                put(2, 1 => H),
                put(2, 1 => Rx(π/4)),
                control(2, 1, 2 => X),
                put(2, 2 => Ry(π/4))
            ),
            obs = ([:Z], [2])
        ),
        (
            name = "Dense gate mix with 2-qubit observable",
            nqubits = 2,
            custom_gates = [
                CliffordGate(:H, [1]),
                PauliRotation(:Y, 1),
                CliffordGate(:CNOT, [1, 2]),
                PauliRotation(:Z, 2)
            ],
            yao_circuit = chain(
                put(2, 1 => H),
                put(2, 1 => Ry(π/4)),
                control(2, 1, 2 => X),
                put(2, 2 => Rz(π/4))
            ),
            obs = ([:X, :Z], [1, 2])
        )
    ]
    for (name, nqubits, custom_gates, yao_circuit, (obs_symbols, obs_qubits)) in circuits
        @testset "$name" begin
            θs = fill(π/4, count(g -> g isa PauliRotation, custom_gates))
            obs = PauliSum(nqubits)
            add!(obs, obs_symbols, obs_qubits, 1.0)
            propagated = propagate(custom_gates, obs, θs)
            custom_val = overlapwithzero(propagated)
            zero_st = zero_state(nqubits)
            evolved_state = Yao.apply(zero_st, yao_circuit)
            yao_obs = build_yao_observable(obs_symbols, obs_qubits, nqubits)
            yao_val = real(expect(yao_obs, evolved_state))
            @test isapprox(custom_val, yao_val; atol=1e-10)
            rev_gates, rev_θs = invert_gates(custom_gates, θs)
            roundtrip_obs = propagate(rev_gates, propagated, rev_θs)
            @test roundtrip_obs == obs
        end
    end
end

@testset "Integration Tests" begin
    # Identity circuit (no gates) on various observables
    for (nqubits, obs_symbols, obs_qubits) in [
            (1, [:X], [1]), 
            (2, [:Z, :X], [1,2]), 
            (3, [:Y, :Z, :X], [1,2,3])
        ]
        @testset "Identity on $nqubits qubit(s), obs=$obs_symbols" begin
            custom_gates = Any[]
            θs = Float64[]
            obs = PauliSum(nqubits)
            add!(obs, obs_symbols, obs_qubits, 1.0)
            propagated = propagate(custom_gates, obs, θs)
            @test propagated == obs
            yao_obs = build_yao_observable(obs_symbols, obs_qubits, nqubits)
            ref_val = real(expect(yao_obs, zero_state(nqubits)))
            custom_val = overlapwithzero(propagated)
            @test isapprox(custom_val, ref_val; atol=1e-10)
        end
    end

    # State overlap functions (overlapwithzero and overlapwithplus)
    @testset "State Overlap Functions" begin
        nqubits = 1
        # Test ⟨0|Z|0⟩ = 1
        obs_z = PauliSum(nqubits)
        add!(obs_z, [:Z], [1], 1.0)
        @test isapprox(overlapwithzero(obs_z), 1.0; atol=1e-10)
        # Test ⟨+|X|+⟩ = 1  via overlapwithplus
        obs_x = PauliSum(nqubits)
        add!(obs_x, [:X], [1], 1.0)
        @test isapprox(overlapwithplus(obs_x), 1.0; atol=1e-10)
        # Test that overlapwithzero(X) = 0 since ⟨0|X|0⟩ = 0
        @test isapprox(overlapwithzero(obs_x), 0.0; atol=1e-10)
        # Test that overlapwithplus(Z) = 0 since ⟨+|Z|+⟩ = 0
        @test isapprox(overlapwithplus(obs_z), 0.0; atol=1e-10)
    end

    # Randomized small circuits with explicit comparison
    @testset "Randomized Circuit Tests" begin
        rng = MersenneTwister(1234)
        for trial in 1:10
            nqubits = rand(rng, 2:3)
            depth = rand(rng, 3:5)
            custom_gates = Any[]
            yao_ops = Any[]
            θs = Float64[]
            for _ in 1:depth
                gate_type = rand(rng, [:H, :X, :Y, :Z, :RX, :RY, :RZ, :CNOT, :SWAP])
                if gate_type == :H
                    q = rand(rng, 1:nqubits)
                    push!(custom_gates, CliffordGate(:H, [q]))
                    push!(yao_ops, put(nqubits, q => H))
                elseif gate_type == :X
                    q = rand(rng, 1:nqubits)
                    push!(custom_gates, CliffordGate(:X, [q]))
                    push!(yao_ops, put(nqubits, q => X))
                elseif gate_type == :Y
                    q = rand(rng, 1:nqubits)
                    push!(custom_gates, CliffordGate(:Y, [q]))
                    push!(yao_ops, put(nqubits, q => Y))
                elseif gate_type == :Z
                    q = rand(rng, 1:nqubits)
                    push!(custom_gates, CliffordGate(:Z, [q]))
                    push!(yao_ops, put(nqubits, q => Z))
                elseif gate_type == :RX
                    q = rand(rng, 1:nqubits)
                    angle = rand(rng) * π/2
                    push!(custom_gates, PauliRotation(:X, q))
                    push!(θs, angle)
                    push!(yao_ops, put(nqubits, q => Rx(angle)))
                elseif gate_type == :RY
                    q = rand(rng, 1:nqubits)
                    angle = rand(rng) * π/2
                    push!(custom_gates, PauliRotation(:Y, q))
                    push!(θs, angle)
                    push!(yao_ops, put(nqubits, q => Ry(angle)))
                elseif gate_type == :RZ
                    q = rand(rng, 1:nqubits)
                    angle = rand(rng) * π/2
                    push!(custom_gates, PauliRotation(:Z, q))
                    push!(θs, angle)
                    push!(yao_ops, put(nqubits, q => Rz(angle)))
                elseif gate_type == :CNOT && nqubits ≥ 2
                    c, t = rand(rng, 1:nqubits), rand(rng, 1:nqubits)
                    while t == c
                        t = rand(rng, 1:nqubits)
                    end
                    push!(custom_gates, CliffordGate(:CNOT, [c, t]))
                    push!(yao_ops, control(nqubits, c, t => X))
                elseif gate_type == :SWAP && nqubits ≥ 2
                    q1, q2 = rand(rng, 1:nqubits), rand(rng, 1:nqubits)
                    while q2 == q1
                        q2 = rand(rng, 1:nqubits)
                    end
                    push!(custom_gates, CliffordGate(:SWAP, [q1, q2]))
                    push!(yao_ops, swap(nqubits, q1, q2))
                end
            end
            if nqubits == 2 && rand(rng) < 0.5
                obs_symbols, obs_qubits = ([:X, :Y], [1, 2])
            else
                q = rand(rng, 1:nqubits)
                sym = rand(rng, [:X, :Y, :Z])
                obs_symbols, obs_qubits = ([sym], [q])
            end
            @testset "Random circuit trial $trial (n=$nqubits, depth=$depth)" begin
                obs = PauliSum(nqubits)
                add!(obs, obs_symbols, obs_qubits, 1.0)
                propagated = propagate(custom_gates, obs, θs)
                custom_val = overlapwithzero(propagated)
                zero_st = zero_state(nqubits)
                evolved_state = apply(zero_st, chain(yao_ops...))
                yao_obs = build_yao_observable(obs_symbols, obs_qubits, nqubits)
                yao_val = real(expect(yao_obs, evolved_state))
                @test isapprox(custom_val, yao_val; atol=1e-10)
                rev_gates, rev_θs = invert_gates(custom_gates, θs)
                roundtrip_obs = propagate(rev_gates, propagated, rev_θs)
                orig_val_zero = overlapwithzero(obs)
                roundtrip_val_zero = overlapwithzero(roundtrip_obs)
                @test isapprox(roundtrip_val_zero, orig_val_zero; atol=1e-10)
            end
        end
    end
end
