### symmetries.jl
##
# This file contains functions to merge Pauli strings by symmetries.
# Currently it supports translational symmetry in 1d and 2d.
##
###


"""
    symmetrymerge(psum::PauliSum, mapfunc::Function)

General symmetry-merging for a PauliSum `psum` using `mapfunc` to transform each 
Pauli string into its representative. 
The expected format for that function is pstr -> mapfunc(pstr).
"""
function symmetrymerge(psum::PauliSum, mapfunc::F) where F<:Function
    merged_psum = PauliPropagation.similar(psum)

    for (pstr, coeff) in psum
        pstr = mapfunc(pstr)
        add!(merged_psum, pstr, coeff)
    end

    return merged_psum
end


"""
    translationmerge(psum::PauliSum)

Shift and merge of a `psum` in a system with 1D translational symmetry.
```
psum = PauliSum(6)
add!(psum, :Z, 3)
add!(psum, :Z, 6)
translationmerge(psum)
>>> PauliSum(nqubits: 6, 1 Pauli term: 
 2.0 * ZIIIII
)
```
"""
translationmerge(psum::PauliSum) = symmetrymerge(
    psum, (pstr -> _translatetolowestinteger(pstr, psum.nqubits))
)


function translationmerge(psum::PauliSum, nx::Integer, ny::Integer)
    if psum.nqubits != nx * ny
        throw(
            ArgumentError("Number of qubits $(psum.nqubits) does not 
                match grid size $(nx) x $(ny)"
            )
        )
    end

    # precompute masks once to accelerate shifting
    main_mask, wrap_mask = _computeshiftleftmasks(paulitype(psum), nx, ny)

    mergefunc(pstr) = _translatetolowestinteger(
        pstr, nx, ny, main_mask, wrap_mask
    )

    return symmetrymerge(psum, mergefunc)
end

function _computeshiftleftmasks(::Type{TT}, nx::Integer, ny::Integer) where TT
    main_mask = zero(TT)
    wrap_mask = zero(TT)

    for col in 1:nx
        for row in 1:ny
            site_index = (row - 1) * nx + col
            bit_index = 2 * (site_index - 1)

            if col == 1
                # first column
                wrap_mask |= (TT(3) << bit_index)
            else
                main_mask |= (TT(3) << bit_index)
            end
        end
    end

    return main_mask, wrap_mask
end

# a function for 1D symmetric merging that does not check for existing terms
# and instead shifts through to find the lowest integer representation
# that is the representative that we merge to
function _translatetolowestinteger(pstr::PauliStringType, nq)
    if pstr == 0
        return pstr
    end

    lowest_pstr = pstr
    for ii in 1:nq
        # shift periodically by one
        pstr = _periodicshiftright(pstr, nq)

        # if the shifted Pauli is lower, break shifting loop
        if pstr < lowest_pstr
            lowest_pstr = pstr
        end
    end

    return lowest_pstr
end

# the same strategy for the 2D case
function _translatetolowestinteger(pstr::PauliStringType, nx, ny, main_mask, wrap_mask)
    if pstr == 0
        return pstr
    end

    lowest_pstr = pstr
    for ii in 1:ny
        for jj in 1:nx
            # shift periodically by one column
            pstr = _periodicshiftleft(pstr, nx, ny, main_mask, wrap_mask)

            # if the shifted Pauli is lower, break shifting loop
            if pstr < lowest_pstr
                lowest_pstr = pstr
            end
        end
        # shift periodically by one row
        if ii < ny
            # shift up one row
            pstr = _periodicshiftup(pstr, nx, ny)
        end
    end

    return lowest_pstr
end

function _periodicshiftright(pstr::PauliStringType, nq)
    first_pauli = getpauli(pstr, 1)
    pstr = PauliPropagation._paulishiftright(pstr)
    pstr = setpauli(pstr, first_pauli, nq)

    return pstr
end



# Shifts a `pstr` left one column in a (`nx`, `ny`) 2D grid of `nq` qubits.
# This function shifts the entire bitstring left one column, 
# and sets the first column of Paulis to the last column of the Paulis.
function _periodicshiftleft(pstr::PauliStringType, nx, ny, main_mask, wrap_mask)
    # main_mask: mask for all bits except the first column
    # wrap_mask: mask for the first column

    shift_size = 2 * nx - 2
    first_col_paulis = pstr & wrap_mask
    main_shift = (pstr & main_mask) >> 2
    pstr = main_shift | (first_col_paulis << shift_size)

    return pstr
end


# Shifts a `pstr` up one row in a (`nx`, `ny`) 2D grid of `nq` qubits on a 
# cylindrical lattice.
# This function shifts the entire bitstring up one row, 
# and sets the first row of Paulis to the last row of the Paulis.
function _periodicshiftup(pstr::PauliStringType, nx, ny)

    n_bits = nx * ny * 2
    shift_size = 2 * nx
    first_row_paulis = pstr & _pauliwindowmask(typeof(pstr), 1, nx)
    pstr = pstr >> shift_size
    pstr = pstr | (first_row_paulis << (n_bits - shift_size))

    return pstr
end

