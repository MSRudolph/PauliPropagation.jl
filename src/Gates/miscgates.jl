### miscgates.jl
##
# A file for gates that don't fit into the other categories.
##
###

using LinearAlgebra

## T gate

struct TGate <: StaticGate
    qind::Int

    @doc """
        TGate(qind::Integer)

    Returns a T gate acting on qubit `qind`.
    It acts on qubit `qind` like a `PauliRotation(:Z, qind)` with angle pi/4.
    """
    TGate(qind::Integer) = (_qinds_check(qind); new(qind))

end

"""
    tomatrix(gate::TGate)

Compute the unitary matrix for a `TGate`.
The returned unitary is returned in Schrödinger picture form. 
"""
function tomatrix(::TGate)
    return _tgate_unitary
end

const _tgate_unitary = [[1 0]; [0 exp(1.0im * pi / 4)]]


## TransferMapGate
# TODO: this should all be made immutable for performance
"""
    TransferMapGate(transfer_map::Vector{Vector{Tuple{PauliStringType,CoeffType}}}, qinds::Vector{Int})
    TransferMapGate(mat::AbstractMatrix, qinds)

A non-parametrized `StaticGate` defined by a transfer map acting on the qubits `qinds`.
Transfer maps can be constructed manually or generated via `totransfermap()`.
"""
struct TransferMapGate{TT,CT} <: StaticGate
    transfer_map::Vector{Vector{Tuple{TT,CT}}}
    qinds::Vector{Int}

    function TransferMapGate(transfer_map::Vector{Vector{Tuple{TT,CT}}}, qinds) where {TT,CT}
        # accept anything that can be converted to a vector of integers
        qinds = vec(collect(qinds))
        nq = length(qinds)

        @assert nq == Int(log(4, length(transfer_map))) "The length of `qinds` `n=$nq` does not match the length of the transfer map `$(length(transfer_map)) ≠ 2^$nq`."

        return new{TT,CT}(transfer_map, qinds)
    end
end


"""
A constructor for `TransferMapGate` that accepts matrix representations in the 0/1 basis or the Pauli basis (a PTM).
"""
function TransferMapGate(mat::AbstractMatrix, qinds)
    # turns number or tuple of numbers into vector of numbers
    qinds = vec(collect(qinds))
    # number of qubits acted on
    nq = length(qinds)

    # infer from the size of the matrix and nq whether it is a matrix in the 0/1 basis or the Pauli basis
    mat_size = size(mat)

    if mat_size != (2^nq, 2^nq) && mat_size != (4^nq, 4^nq)
        throw(ArgumentError("The matrix must be square and have size (2^$nq x 2^$nq) or (4^$nq x 4^$nq) " *
                            "given the passed qinds=$qinds."))
    end

    if mat_size == (2^nq, 2^nq)
        # the matrix is assumed to be in the 0/1 basis
        # transform it into a PTM
        mat = calculateptm(mat)
    end

    ptmap = totransfermap(mat)
    return TransferMapGate(ptmap, qinds)
end

# Reconstruct 4×4 PTM from 1-qubit transfer map (sparse column-wise representation).
function _transfermap_to_ptm_1qubit(transfer_map)
    ptm = zeros(4, 4)
    for j in 1:4
        for (pstr, coeff) in transfer_map[j]
            ptm[Int(pstr) + 1, j] = coeff
        end
    end
    return ptm
end

# Reconstruct 2×2 unitary from 1-qubit PTM (up to global phase).
# PTM[i,j] = tr(P_i U P_j U†). The 3×3 block R = ptm[2:4, 2:4] is the Bloch-sphere rotation.
function _ptm_to_unitary_1qubit(ptm)
    R = ptm[2:4, 2:4]
    # Rotation angle: trace(R) = 1 + 2cos(theta)
    cos_theta = (tr(R) - 1) / 2
    cos_theta = clamp(cos_theta, -1.0, 1.0)
    theta = acos(cos_theta)

    if abs(sin(theta)) < 1e-12
        # theta ≈ 0 or pi: identity or 180° rotation; axis arbitrary for theta=0
        if theta < 1e-10
            return Matrix{ComplexF64}(I, 2, 2)
        end
        # theta ≈ pi: axis from eigenvector for eigenvalue 1
        F = eigen(R)
        idx = findfirst(x -> isapprox(x, 1.0), F.values)
        n = vec(real.(F.vectors[:, idx]))
        n = n / norm(n)
    else
        n_x = (R[2, 3] - R[3, 2]) / (2 * sin(theta))
        n_y = (R[3, 1] - R[1, 3]) / (2 * sin(theta))
        n_z = (R[1, 2] - R[2, 1]) / (2 * sin(theta))
        n = [n_x, n_y, n_z]
    end

    # U = cos(theta/2) I - i sin(theta/2) (n·σ)
    c = cos(theta / 2)
    s = sin(theta / 2)
    nx, ny, nz = n[1], n[2], n[3]
    return [
        c - 1.0im*s*nz    -1.0im*s*(nx - 1.0im*ny);
        -1.0im*s*(nx + 1.0im*ny)   c + 1.0im*s*nz
    ]
end

"""
    tomatrix(gate::TransferMapGate)

Return the unitary matrix for a `TransferMapGate` in the computational basis.
Supported only for single-qubit gates: the unitary is reconstructed from the
transfer map (up to global phase). Multi-qubit gates do not support `tomatrix`.
"""
function tomatrix(gate::TransferMapGate)
    nq = length(gate.qinds)
    if nq != 1
        throw(ArgumentError(
            "tomatrix(TransferMapGate) is only supported for single-qubit gates, got $nq qubits."
        ))
    end
    ptm = _transfermap_to_ptm_1qubit(gate.transfer_map)
    return _ptm_to_unitary_1qubit(ptm)
~end