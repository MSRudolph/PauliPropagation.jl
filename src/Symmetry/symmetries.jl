### symmetries.jl
##
# This file contains functions to merge Pauli strings by symmetries.
# Currently it supports translational symmetry in 1d and 2d.
##
###

#TODO(YT): remove the commented code below if Manuel is happy with the short hand
# """
#     translationmerging(psum::PauliSum)
# Shift and merge of a `psum` in a system with translational symmetry.
# ```
# psum = PauliSum(6)
# add!(psum, :Z, 3)
# add!(psum, :Z, 6)
# translationmerging(psum)
# >>> PauliSum(nqubits: 6, 1 Pauli term: 
#  2.0 * ZIIIII
# )
# ```
# """
# function translationmerging(psum::PauliSum)
#     shifted_psum = PauliPropagation.similar(psum)

#     for (pstr, coeff) in psum
#         pstr = _translatetolowestinteger(pstr, shifted_psum.nqubits)
#         add!(shifted_psum, pstr, coeff)
#     end

#     return shifted_psum
# end


function translationmerging(psum::PauliSum, nx, ny)
    shifted_psum = PauliPropagation.similar(psum)

    if shifted_psum.nqubits != nx * ny
        throw(
            ArgumentError("Number of qubits $(shifted_psum.nqubits) does not \n
                match grid size $(nx) x $(ny)"
            )
        )
    end

    for (pstr, coeff) in psum
        pstr = _translatetolowestinteger(pstr, nx, ny)
        add!(shifted_psum, pstr, coeff)
    end

    return shifted_psum
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
        pstr = _periodictshiftright(pstr, nq)

        # if the shifted Pauli is lower, break shifting loop
        if pstr < lowest_pstr
            lowest_pstr = pstr
        end
    end

    return lowest_pstr
end

# the same strategy for the 2D case
function _translatetolowestinteger(pstr::PauliStringType, nx, ny)
    if pstr == 0
        return pstr
    end

    lowest_pstr = pstr
    for ii in 1:ny
        for jj in 1:nx
            # shift periodically by one column
            pstr = _periodicshiftleft(pstr, nx, ny)

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

function _periodictshiftright(pstr::PauliStringType, nq)
    first_pauli = getpauli(pstr, 1)
    pstr = PauliPropagation._paulishiftright(pstr)
    pstr = setpauli(pstr, first_pauli, nq)

    return pstr
end



# Shifts a `pstr` left one column in a (`nx`, `ny`) 2D grid of `nq` qubits.
# This function shifts the entire bitstring left one column, 
# and sets the first column of Paulis to the last column of the Paulis.
function _periodicshiftleft(pstr::PauliStringType, nx, ny)
    # TODO: optimize setpauli for non-contiguous indices. 
    # But Manuel has a mask function that can be used for this purpose!

    nq = nx * ny

    col_paulis = getpauli(pstr, 1:nx:nq) # get the first column of Paulis
    pstr = pstr >> 2
    pstr = setpauli(pstr, col_paulis, nx:nx:nq)

    return pstr
end


# Shifts a `pstr` up one row in a (`nx`, `ny`) 2D grid of `nq` qubits on a 
# cylindrical lattice.
# This function shifts the entire bitstring up one row, 
# and sets the first row of Paulis to the last row of the Paulis.
function _periodicshiftup(pstr::PauliStringType, nx, ny)

    nq = nx * ny

    row_paulis = getpauli(pstr, 1:nx) # get the first row of Paulis
    pstr = pstr >> (2 * nx)  # shift by nx many Paulis
    pstr = setpauli(pstr, row_paulis, (nq-nx+1):nq)
    # set the last row of Paulis to the first row of the bitstring

    return pstr
end

"""
    greedytranslationmerging(pstr::PauliStringType)

Shift a `pstr` to the right until the first non-I Pauli, used for a system with
translational symmetry in 1d.
"""
function greedytranslationmerging(pstr::PauliStringType, nq)
    if pstr == 0
        return pstr
    end

    # shift until the first non-zero Pauli
    while getpauli(pstr, 1) == 0
        pstr = _periodictshiftright(pstr, nq)
    end

    return pstr
end

"""
    generalmerging(psum::PauliSum, mergefunc::Function)
General merging of a `psum` using a user defined `mergefunc`.
"""
function generalmerging(psum::PauliSum, mergefunc::Function)
    merged_psum = PauliPropagation.similar(psum)

    for (pstr, coeff) in psum
        pstr = mergefunc(pstr, psum.nqubits)
        add!(merged_psum, pstr, coeff)
    end

    return merged_psum
end

"""
    greedytranslationmerging(psum::PauliSum)

Greedy translation and merge of a `psum` in 1d system.
"""
greedytranslationmerging(psum::PauliSum) = generalmerging(
    psum, (pstr, nq) -> greedytranslationmerging(pstr, nq)
)


"""
    translationmerging(psum::PauliSum)

Shift and merge of a `psum` in a system with translational symmetry.
```
psum = PauliSum(6)
add!(psum, :Z, 3)
add!(psum, :Z, 6)
translationmerging(psum)
>>> PauliSum(nqubits: 6, 1 Pauli term: 
 2.0 * ZIIIII
)
```
"""
translationmerging(psum::PauliSum) = generalmerging(
    psum, (pstr, nq) -> _translatetolowestinteger(pstr, nq)
)
