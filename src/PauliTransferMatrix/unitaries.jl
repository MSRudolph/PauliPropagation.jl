###
##
# This file contains implementations of certain unitary gates.
# Support gates are local unitaries, such as one-qubit, two-qubit Pauli gates.
# Note our convetion is to define the conjugate transpose of the unitary matrix.
# In the propagation, the circuit is ran in the reverse order, so defining the
# conjugate transpose of the unitary matrix is more convenient.
##
###

# module unitaries

using LinearAlgebra

struct U1Unitary <: ParametrizedGate
    nparams::Int
    qinds::Vector{Int}
end

struct TUnitary <: StaticGate
    qinds::Vector{Int}
end

const Idmat = [1 0; 0 1]
const Xmat = [0 1; 1 0]
const Ymat = [0 -1im; 1im 0]
const Zmat = [1 0; 0 -1]

const pauli_basis = [Idmat / sqrt(2), Xmat / sqrt(2), Ymat / sqrt(2), Zmat / sqrt(2)]

"""
    function constant_pauli_basis_n_qubits(n::Int)

Compute the Pauli basis for `n` qubits.

Arguments
- `n::Int`: The number of qubits.

Returns
- `basis::Vector{Array{ComplexF64}}`: The Pauli basis for `n` qubits.
"""
function constant_pauli_basis_n_qubits(n::Int)
    if n == 1
        return pauli_basis
    else
        # Compute the Kronecker product iteratively for n qubits
        basis = pauli_basis
        for _ in 2:n
            basis = [kron(p1, p2) for p1 in basis, p2 in pauli_basis]
        end
        return basis
    end
end

"""
    function get_unitary_dagmat(gate::TUnitary, params)

Compute the conjugate transpose of the unitary matrix for the T gate.

Arguments
- `gate::TUnitary`: The T gate.
- `params::Vector{Float64}`: The parameters for the gate.

Returns
- `U::Array{ComplexF64}`: The conjugate transpose of the unitary matrix.
"""
function get_unitary_dagmat(gate::TUnitary, params=[])

    U = [[1 0]; [0 exp(1.0im * pi / 4)]]

    return U'
end

"""
    function get_unitary_dagmat(gate::PauliRotation, params)

Compute the conjugate transpose of the unitary matrix for the PauliRotation gate.

Arguments
- `gate::PauliRotation`: The PauliRotation gate.
- `params::Vector{Float64}`: A length-1 parameter for the gate.

Returns 
- `U::Array{ComplexF64}`: The conjugate transpose of the unitary matrix.
"""    
function get_unitary_dagmat(gate::PauliRotation, params)

    theta = params[1]
    nqubits = length(gate.qinds)
    pauli_basis_vec = constant_pauli_basis_n_qubits(nqubits)

    # The indices of the pauli matrix are sorted in ascending order
    sorted_indices = sortperm(gate.qinds)
    pauli = gate.symbols[sorted_indices]

    # These pauli matrices are normalized by sqrt(2)^nqubits
    pauli_mat = sqrt(2)^nqubits * pauli_basis_vec[symboltoint(pauli) + 1]

    id = Matrix(1.0I, 2^nqubits, 2^nqubits)

    U = cos(theta/2) * id - 1.0im * sin(theta/2) * pauli_mat

    return U'
end

"""
    function get_unitary_dagmat(gate::U1Unitary, params)

Compute the conjugate transpose of the unitary matrix for the U1 gate.

The 4x4 unitary is of the form:
U = [[exp(iφ1) 0 0 0];
     [0 exp(-0.5im * (α + β - 2δ)) * cos(θ/2) exp(0.5im * (α - β + 2δ)) * sin(θ/2) 0];
     [0 - exp( - 0.5im * (α - β - 2δ)) * sin(θ/2) exp(0.5im * (α + β + 2δ)) * cos(θ/2) 0];
     [0 0 0 exp(iφ2)]]

Arguments
- `gate::U1Unitary`: The U1 gate.
- `params::Vector{Float64}`: The parameters for the gate.
    Order of the parameters are [θ, α, β, δ, φ1, φ2].

Returns
- `U::Array{ComplexF64}`: The conjugate transpose of the unitary matrix.
"""
function get_unitary_dagmat(gate::U1Unitary, params)

    nqubits = length(gate.qinds)
    all_params = zeros(Float64, gate.nparams)
    all_params[1:length(params)] = params
    θ, α, β, δ, φ1, φ2 = all_params

    U = zeros(ComplexF64, 2^nqubits, 2^nqubits)

    # Add checks for the number of qubits
    U[1, 1] = exp(1.0im * φ1)
    U[4, 4] = exp(1.0im * φ2)
    U[2, 2] = exp(-0.5im * (α + β - 2δ)) * cos(θ/2)
    U[2, 3] = exp(0.5im * (α - β + 2δ)) * sin(θ/2)
    U[3, 2] = - exp( - 0.5im * (α - β - 2δ)) * sin(θ/2)
    U[3, 3] = exp(0.5im * (α + β + 2δ)) * cos(θ/2)

    return U'
end
