using SparseArrays


"""
    function calc_ptm_dagmap(udag; tol=1e-15, symmetry=nothing)

Arguments
- `udag::Array{}`: The unitary matrix for which the PTM is calculated.
- `tol::Float64=1e-15`: The tolerance for dropping small values in the PTM.
- `symmetry::Array{}`: The symmetry matrix to apply to the PTM.
    If no symmetry is provided, the PTM is computed in the Pauli basis.

Returns
- `ptm::SparseMatrixCSC`: The PTM of the unitary matrix `udag`.
"""
function calc_ptm_dagmap(udag; tol=1e-15, symmetry=nothing)

    nqubits = Int(log(2, size(udag)[1]))

    # Write ptm as a sprase matrix
    ptm = spzeros(Float64, 4^nqubits, 4^nqubits)

    # The pauli basis vector is defined to be consistent with index of the pstr
    pauli_basis_vec = constant_pauli_basis_n_qubits(nqubits)

    # We are taking udag as the actual unitary, and the PTM is defined by
    # evolving P_j and take overlap with P_i
    # PTM_{ij} = Tr(udag * P_j * udag^{\dagger} * P_i)
    # in the Pauli basis, this is always real.
    for i in 1:4^nqubits
        for j in 1:4^nqubits
            ptm[i, j] = real(
                tr(udag * pauli_basis_vec[j] * udag' * pauli_basis_vec[i])
            )
        end
    end

    if symmetry != nothing
        ptm = sparse(inv(symmetry) * ptm * symmetry)
    end

    droptol!(ptm, tol)  # Here tol is relevant to the input pstr.

    return ptm  # return the ptm as a sparse matrix
end

"""
    function get_ptm_sparse(ptm::SparseMatrixCSC, type::DataType)

Arguments
- `ptm::SparseMatrixCSC`: The PTM in sparse matrix format.
- `type::Function`: The type to apply to the indices of the PTM.

Returns
- `ptm_map::Vector{Tuple}`: The PTM in a sparse format.
"""
function get_ptm_sparse(ptm::SparseMatrixCSC, type)
    ptm_map = Vector{Tuple}()
    for j in 1:size(ptm, 2) # loop over columns
        rows, vals = findnz(view(ptm, :, j))
        rows = [i - 1 for i in rows] # convert to 0-based indexing used for PP
        # Create a tuple of the form (p1, c1, p2, c2, ...)
        column_tuple = Tuple(Iterators.flatten(zip(type.(rows), vals)))
        push!(ptm_map, column_tuple)
    end
    return ptm_map
end

"""
    function get_ptm_sparse(ptm::Array, type::DataType)

Arguments
- `ptm::Array`: The PTM in matrix format.
- `type::Function`: The type to apply to the indices of the PTM.

Returns
- `ptm_map::Vector{Tuple}`: The PTM in a sparse format.
"""
function get_ptm_sparse(ptm, type)

    ptm_map = Vector{Tuple}()
    rows = collect(range(0, size(ptm, 1) - 1)) # convert to 0-based indexing used for PP

    for j in 1:size(ptm, 2) # loop over columns

        # find non-zero elements
        rows = findall(!iszero, ptm[:, j])
        vals = ptm[rows, j]
        rows_0based = [i - 1 for i in rows] # convert to 0-based indexing used for PP

        push!(ptm_map, Tuple(Iterators.flatten(zip(type.(rows_0based), vals))))

    end
    return ptm_map
end

