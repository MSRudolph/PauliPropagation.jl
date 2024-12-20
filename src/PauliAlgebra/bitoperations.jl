"""
    getinttype(nqubits::Integer)

Function to return the smallest integer type that can hold nqubits for memory and speed.
"""
function getinttype(nqubits::Integer)
    # TODO: This function is type unstable

    # we need 2 bits per qubit
    nbits = 2 * nqubits

    # select correct UInt type
    if nbits <= 8
        inttype = UInt8
    elseif nbits <= 16
        inttype = UInt16
    elseif nbits <= 32
        inttype = UInt32
    elseif nbits <= 64
        inttype = UInt64
    elseif nbits <= 128
        inttype = UInt128
    elseif nbits <= 256    # we have a custom hash function for these
        inttype = UInt256
    elseif nbits <= 1024
        inttype = BigInt  # TODO: get larger Integer types that aren't BigInt
    else
        throw("The maximum number of qubits supported is currently 1024.")
    end

    return inttype
end

import Base: hash, hash_uint
"""
Re-defining the hash of `BitIntegers.UInt256` via type piracy.
"""
Base.hash(x::UInt256) = hash_uint(x)

"""
    hash_uint(x::UInt256)

Custom hash-function for `BitIntegers.UInt256`. 
It appears faster in practice than the default hash for this type, but it is not fully tested.
It certainly allocates no memory. 
"""
function hash_uint(x::UInt256)
    mask::UInt64 = 0xFFFFFFFFFFFFFFFF
    # hash each 64-bit segment using Julia's inbuilt hashes
    a::UInt64 = hash_uint(UInt64((x >> 0) & mask))
    b::UInt64 = hash_uint(UInt64((x >> 64) & mask))
    c::UInt64 = hash_uint(UInt64((x >> 128) & mask))
    d::UInt64 = hash_uint(UInt64((x >> 192) & mask))
    # mess it up a little using this ChatGPT hallucination
    h::UInt64 = 0
    h = a ⊻ ((b << 21) | (b >> (64 - 21)))
    h = h ⊻ ((c << 42) | (c >> (64 - 42)))
    h = h ⊻ ((d << 63) | (d >> (64 - 63)))
    h = (h * 0x9e3779b185ebca87) ⊻ (h >> 32)
    h = (h * 0xc2b2ae3d27d4eb4f) ⊻ (h >> 29)
    h = (h * 0x165667b19e3779f9) ⊻ (h >> 32)
    return h
end

"""
    _countbitweight(pstr::PauliStringType)

This function counts the number of 00 bit pairs in the integer Pauli string.
"""
function _countbitweight(pstr::PauliStringType)
    # get our super bit mask looking like ....1010101.
    mask = alternatingmask(pstr)

    # m1 carries the 1's of the Paulis on odd bits 
    m1 = pstr & mask

    # m2 carries the 1's of the pauliP on even bits
    m2 = pstr & (mask << 1)

    # OR between m1 and left-shifted m2 to get 1's where either m1 or m2 arre 1
    res = m1 | (m2 >> 1)

    # count 1's to get the number of non-identity Paulis
    return count_ones(res)
end

"""
     _countbitxy(pstr::PauliStringType)

This function counts the number of 01 (X) or 10 (Y) bit pairs in the integer Pauli string.
"""
function _countbitxy(pstr::PauliStringType)
    # we use that 01 and 10 have exactly one 1 and one 0

    # get our super bit mask looking like ....1010101.
    mask = alternatingmask(pstr)

    # XOR to put 1's where the bits are different
    op = pstr ⊻ (pstr >> 1)

    # AND with the mask to extract the 1's
    op = op & mask

    # count 1's to get the number of X or Y Paulis
    return count_ones(op)
end

"""
    _countbityz(pstr::PauliStringType)

This function counts the number of 10 (Y) or 11 (Z) bit pairs in the integer Pauli string.
"""
function _countbityz(pstr::PauliStringType)
    # we use that both have a 1 on the left bit

    # get our super bit mask looking like ....1010101.
    mask = alternatingmask(pstr)

    # AND with the shifted mask to extract the 1's on the left bit
    op = pstr & (mask << 1)

    # count 1's to get the number of Y or Z Paulis
    return count_ones(op)
end

"""
    _bitcommutes(pstr1::PauliStringType, pstr2::PauliStringType)

This function checks if two integer Pauli strings commute.
"""
function _bitcommutes(pstr1::PauliStringType, pstr2::PauliStringType)

    mask0 = alternatingmask(pstr1)
    mask1 = mask0 << 1

    # obtain the left (then right) bits of each Pauli pair, in-place
    aBits0 = mask0 & pstr1
    aBits1 = mask1 & pstr1
    bBits0 = mask0 & pstr2
    bBits1 = mask1 & pstr2

    # shift left bits to align with right bits
    aBits1 = aBits1 >> 1
    bBits1 = bBits1 >> 1

    # sets '10' at every Pauli index where individual pairs don't commute
    flags = (aBits0 & bBits1) ⊻ (aBits1 & bBits0)

    # strings commute if parity of non-commuting pairs is even
    return (count_ones(flags) % 2) == 0
end

"""
    _bitpaulimultiply(pstr1::PauliStringType, pstr2::PauliStringType)

XOR between two Pauli different non-identity strings gives the third one. Ignores signs or any coefficient.
"""
_bitpaulimultiply(pstr1::PauliStringType, pstr2::PauliStringType) = pstr1 ⊻ pstr2

"""
    _paulishiftright(pstr::PauliStringType)

Shift to the right and truncate the first encoded Pauli string. Just a utility function.
"""
_paulishiftright(pstr::PauliStringType) = pstr >> 2


"""
    _getpaulibits(pstr::PauliStringType, index::Integer)

This function extracts the Pauli at position `index` from the integer Pauli string.
"""
function _getpaulibits(pstr::PauliStringType, index::Integer)
    # we need to shift the integer by 2 * (index - 1), then the first two bits are target Pauli
    bitindex = 2 * (index - 1)

    # shift to the right
    shifted_pstr = (pstr >> bitindex)

    # AND with 3 (00000011) to get the first two bits
    return shifted_pstr & typeof(pstr)(3)
end

"""
    _getbit(pauli::Integer, bitindex::Integer)

Gets the bit at index `bitindex` in the integer Pauli string.
"""
function _getbit(pauli::Integer, bitindex::Integer)
    # return integer with ...000[bit].

    # shift by bitindex
    shifted_pauli = (pauli >> bitindex)

    # AND with 1 to get first bit
    return shifted_pauli & typeof(pauli)(1)
end

"""
    _setpaulibits(pstr::PauliStringType, pauli::PauliType, index::Integer)

This function sets the Pauli at position `index` in the integer Pauli string to `target_pauli`.
"""
function _setpaulibits(pstr::PauliStringType, target_pauli::PauliType, index::Integer)
    # we need to shift the integer by 2 * (index - 1), then the first two bits are target Pauli
    bitindex = 2 * (index - 1)

    # read bits of the pauli
    b1 = _getbit(target_pauli, 0)
    b2 = _getbit(target_pauli, 1)

    # insert them into the pstr
    pstr = _setbit(pstr, b1, bitindex)
    pstr = _setbit(pstr, b2, bitindex + 1)
    return pstr
end

"""
    _setpaulibits(
        pstr::PauliStringType, 
        target_pauli::PauliType, 
        indices::Vector{Int}
    )

Sets multiple Paulis at `indices` in `pstr` to `target_pauli`.
"""
function _setpaulibits(
    pstr::PauliStringType,
    target_pauli::PauliType,
    indices::Vector{Int}
)
    # set multiple Paulis at once
    for (ii, index) in enumerate(indices)
        bitindex = 2 * (index - 1)

        # read bits of the pauli
        pauli = _getpaulibits(target_pauli, ii)
        b1 = _getbit(pauli, 0)
        b2 = _getbit(pauli, 1)

        # insert them into the pstr
        pstr = _setbit(pstr, b1, bitindex)
        pstr = _setbit(pstr, b2, bitindex + 1)
    end
    return pstr
end

"""
    _setbit(pstr::PauliStringType, target_bit::Integer, bitindex::Integer)

Sets a bit at index `bitindex` in the integer Pauli string to the value of `target_bit`.
"""
function _setbit(pstr::PauliStringType, target_bit::Integer, bitindex::Integer)
    # set bit at bitindex to bit

    if target_bit == true  # set to one
        pstr = _setbittoone(pstr, bitindex)
    else
        pstr = _setbittozero(pstr, bitindex)
    end
    return pstr
end

"""
    _setbittoone(pstr::Integer, bitindex::Integer)

Sets a bit at index `bitindex` in the integer Pauli string to 1.
"""
function _setbittoone(pstr::Integer, bitindex::Integer)
    # set bit at bitindex to 1

    # shift ...00100...  to bitindex
    shifted_onebit = (typeof(pstr)(1) << bitindex)

    # OR with pauli string to make sure that that bit is 1
    return pstr | shifted_onebit
end


"""
    _setbittozero(pstr::Integer, bitindex::Integer)

Sets a bit at index `bitindex` in the integer Pauli string to 0.
"""
function _setbittozero(pstr::Integer, bitindex::Integer)
    # set bit at bitindex to 0

    # flip all bits
    pstr = ~pstr

    # set target bit to one
    pstr = _setbittoone(pstr, bitindex)

    # flip all bits back, only the target bit is 0
    pstr = ~pstr
    return pstr
end


@generated function alternatingmask(pstr::T) where {T<:PauliStringType}
    # define our super bit mask looking like ....1010101.

    # length is the number of bits in the integer
    n_bits = min(bitsize(pstr), 2_048)  # for max 1024 qubits.
    mask = pstr(0)
    for ii in 0:(n_bits-1)
        if ii % 2 == 0
            mask = _setbittoone(mask, ii)
        end
    end
    return mask
end