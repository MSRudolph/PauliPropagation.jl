using Test
using LinearAlgebra
using Yao
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
