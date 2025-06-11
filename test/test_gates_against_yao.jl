using Test
using LinearAlgebra
using Random
using Yao: X, Y, Z, H, Rx, Rz, Ry, chain, put, control, zero_state, expect, apply, rot, mat, matblock, swap, SWAP
using PauliPropagation

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
const single_obs = [:X, :Y, :Z]
const two_obs = [(:X, :X), (:X, :Y), (:X, :Z),
                  (:Y, :X), (:Y, :Y), (:Y, :Z),
                  (:Z, :X), (:Z, :Y), (:Z, :Z)]

@testset "Observable Propagation Check via Yao.jl" begin
    for gate_symbol in all_gates
        nqubits = gate_symbol in (:CNOT, :CZ, :SWAP, :ZZpihalf) ? 2 : 1
        yao_gate = clifford_to_yao(CliffordGate(gate_symbol, 1:nqubits))
        for obs_symbol in single_obs
            @testset "Gate $gate_symbol propagates 1-local $obs_symbol" begin
                qc = [CliffordGate(gate_symbol, 1:nqubits)]
                obs = PauliSum(nqubits)
                add!(obs, [obs_symbol], [1], 1)
                propagated = propagate(qc, obs)
                exp_val = overlapwithzero(propagated)

                state = zero_state(nqubits)
                evolved = apply(state, yao_gate)
                yao_obs = put(nqubits, 1 => (obs_symbol == :X ? X : obs_symbol == :Y ? Y : Z))
                ref_val = real(expect(yao_obs, evolved))
                @test isapprox(exp_val, ref_val; atol=1e-10)
            end
        end
        if nqubits == 2
            for (o1, o2) in two_obs
                @testset "Gate $gate_symbol propagates 2-local $o1⊗$o2" begin
                    qc = [CliffordGate(gate_symbol, 1:2)]
                    obs = PauliSum(2)
                    add!(obs, [o1, o2], [1, 2], 1)
                    propagated = propagate(qc, obs)
                    exp_val = overlapwithzero(propagated)
                    state = zero_state(2)
                    evolved = apply(state, yao_gate)
                    p1 = o1 == :X ? X : o1 == :Y ? Y : Z
                    p2 = o2 == :X ? X : o2 == :Y ? Y : Z
                    yao_obs = chain(put(2, 1 => p1), put(2, 2 => p2))
                    ref_val = real(expect(yao_obs, evolved))
                    @test isapprox(exp_val, ref_val; atol=1e-10)
                end
            end
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

function register_inverse!(symbol::Symbol)
    inv_symbol = Symbol(string(symbol), "_inv")
    if !haskey(clifford_map, inv_symbol)
        original_map = clifford_map[symbol]
        clifford_map[inv_symbol] = transposecliffordmap(original_map)
    end
    return inv_symbol
end

function invert_gates(gates::Vector{<:Any}, θs::Vector{Float64})
    rev_list = Any[]
    rev_θs   = Float64[]
    θ_index = length(θs)
    for gate in reverse(gates)
        if gate isa PauliRotation
            push!(rev_list, gate)
            push!(rev_θs, -θs[θ_index])
            θ_index -= 1
        elseif gate isa CliffordGate
            inv_sym = register_inverse!(gate.symbol)
            push!(rev_list, CliffordGate(inv_sym, gate.qinds))
        else
            throw(ArgumentError("Cannot invert object of type $(typeof(gate))"))
        end
    end
    return rev_list, rev_θs
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
        )
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
            yao_obs = build_yao_observable(obs_symbols, obs_qubits, circuit.nqubits)
            yao_val = real(expect(yao_obs, evolved_state))
            @test isapprox(custom_val, yao_val; atol=1e-10)
            rev_gates, rev_θs = invert_gates(circuit.custom_gates, θs)
            roundtrip_obs = propagate(rev_gates, propagated, rev_θs)
            @test roundtrip_obs == obs
        end
    end
end

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
