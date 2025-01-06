using LinearAlgebra
using Random
using Test

tol=1e-10

@testset "Pauli basis" begin
    """Test recursive generation of Pauli basis for 1 and 2 qubits"""
    nqubits = 1
    @test constant_pauli_basis_n_qubits(nqubits) == pauli_basis

    nqubits = 2
    pauli_basis_two_qubit = [
        kron(p1, p2) for p1 in pauli_basis, p2 in pauli_basis
    ]
    @test constant_pauli_basis_n_qubits(nqubits) == pauli_basis_two_qubit
    
end

@testset "U1 unitary two-qubit" begin
    
    gate = U1Unitary(6, [0, 1])
    params = [0.1, 0.2, 0.3, 0.4, 0., 0.]
    expected_U = [
        [1 0 0 0 ]; 
        [0 0.9875353715596337 + 0.1492513737209447im 0.04694906782365546 + 0.01713774756136047im 0 ]; 
        [0 -0.045003598147776255 - 0.021739216056256186im 0.7950889010970834 + 0.6044300803163632im 0 ]; 
        [0 0 0 1]
    ]  # this is the expected unitary matrix for the given parameters (QuEst)
    U = get_unitary_dagmat(gate, params)'
    
    @test LinearAlgebra.norm(U - expected_U) < tol
end

@testset "PauliRotation two-qubit" begin
    # test Pauli rotation daggered unitary
    pauligate = PauliRotation([:I, :X], [1, 3])
    theta = Random.rand()
    Rxdag = get_unitary_dagmat(pauligate, [theta])
    
    expected_Rxdag = [
        [cos(theta/2) 1im * sin(theta/2) 0 0];
        [1im * sin(theta/2) cos(theta/2) 0 0];
        [0 0 cos(theta/2) 1im * sin(theta/2)];
        [0 0 1im * sin(theta/2) cos(theta/2)]
    ]
    @test LinearAlgebra.norm(Rxdag - expected_Rxdag) < tol
end

@testset "T gate" begin
    # test T unitary
    tgate = TUnitary([1])
    UT = get_unitary_dagmat(tgate, [])
    expected_UT = [
        [1 0];
        [0 exp(-1.0im * pi / 4)]
    ]
    @test LinearAlgebra.norm(UT - expected_UT) < tol
end


