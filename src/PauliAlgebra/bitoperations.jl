# testing type aliases here
const PauliStringType = Union{UInt8,UInt16,UInt32,UInt64,UInt128,UInt256,BigInt,Int} # to be maintained when we adapt getinttype()
const SinglePauliType = PauliStringType


function getinttype(nqubits)
    # get the smallest integer type that can hold nq qubits for memory and speed

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
    elseif nbits <= 256    # These are super slow to hash # TODO: custom hash functions
        inttype = UInt256
        # elseif nbits <= 512
        #     inttype = UInt512
        # elseif nbits <= 1024
        #     inttype = UInt1024
    else
        inttype = BigInt  # TODO: get larger Integer types that aren't BigInt
    end

    return inttype
end

function countbitweight(paulistring_int::PauliStringType; kwargs...)
    # this function counts the number of 00 bit pairs in the integer representation of the Pauli string

    # get our super bit mask looking like ....1010101.
    mask = alternatingmask(paulistring_int)

    # m1 carries the 1's of oper on odd bits 
    m1 = paulistring_int & mask

    # m2 carries the 1's of oper on even bits
    m2 = paulistring_int & (mask << 1)

    # OR between m1 and left-shifted m2 to get 1's where either m1 or m2 arre 1
    res = m1 | (m2 >> 1)

    # count 1's to get the number of non-identity Pauli operators
    return count_ones(res)
end


function countbitxy(paulistring_int::PauliStringType; kwargs...)
    # this function counts the number of 01 (X) or 10 (Y) bit pairs in the integer representation of the Pauli string
    # we use that 01 and 10 have exactly one 1 andd one 0

    # get our super bit mask looking like ....1010101.
    mask = alternatingmask(paulistring_int)

    # XOR to put 1's where the bits are different
    op = paulistring_int ⊻ (paulistring_int >> 1)

    # AND with the mask to extract the 1's
    op = op & mask

    # count 1's to get the number of X or Y operators
    return count_ones(op)
end

function countbityz(paulistring_int::PauliStringType; kwargs...)
    # this function counts the number of 10 (Y) or 11 (Z) bit pairs in the integer representation of the Pauli string
    # we use that both have a 1 on the left bit

    # get our super bit mask looking like ....1010101.
    mask = alternatingmask(paulistring_int)

    # AND with the shifted mask to extract the 1's on the left bit
    op = paulistring_int & (mask << 1)

    # count 1's to get the number of Y or Z operators
    return count_ones(op)
end

function bitcommutes(paulistring1_int::PauliStringType, paulistring2_int::PauliStringType)
    # TODO: document this nicely

    mask0 = alternatingmask(paulistring1_int)
    mask1 = mask0 << 1

    # obtain the left (then right) bits of each Pauli pair, in-place
    aBits0 = mask0 & paulistring1_int
    aBits1 = mask1 & paulistring1_int
    bBits0 = mask0 & paulistring2_int
    bBits1 = mask1 & paulistring2_int

    # shift left bits to align with right bits
    aBits1 = aBits1 >> 1
    bBits1 = bBits1 >> 1

    # sets '10' at every Pauli index where individual pairs don't commute
    flags = (aBits0 & bBits1) ⊻ (aBits1 & bBits0)

    # strings commute if parity of non-commuting pairs is even
    return (count_ones(flags) % 2) == 0
end

# 'xor between two Pauli different non-identity strings gives the third one. Ignores signs or any coefficient.
bitpaulimultiply(paulistring1_int::PauliStringType, paulistring2_int::PauliStringType) = paulistring1_int ⊻ paulistring2_int

# shift to the right and truncate the first encoded Pauli string. Just a utility function.
_paulishiftright(paulistring_int::PauliStringType) = paulistring_int >> 2

function getsinglepaulibits(paulistring_int::PauliStringType, index::Integer)
    # this function extracts the Pauli operator at position index from the integer representation of the Pauli string

    # we need to shift the integer by 2 * (index - 1), then the first two bits are target Pauli
    bitindex = 2 * (index - 1)

    # shift to the right
    shifted_oper = (paulistring_int >> bitindex)

    # AND with 3 (00000011) to get the first two bits
    return shifted_oper & typeof(paulistring_int)(3)
end


function setsinglepaulibits(paulistring_int::PauliStringType, index::Integer, pauli_int::SinglePauliType)
    # this function sets the Pauli operator at position index in the integer representation of the Pauli string

    # we need to shift the integer by 2 * (index - 1), then the first two bits are target Pauli
    bitindex = 2 * (index - 1)

    # read bits of the pauli_int
    b1 = _readbit(pauli_int, 0)
    b2 = _readbit(pauli_int, 1)

    # insert them into the paulistring_int
    paulistring_int = _setbit(paulistring_int, bitindex, b1)
    paulistring_int = _setbit(paulistring_int, bitindex + 1, b2)
    return paulistring_int
end

function _setbit(paulistring_int::PauliStringType, bitindex::Integer, target_bit::Integer)
    # set bit at bitindex to bit

    if target_bit == true  # set to one
        paulistring_int = _setbittoone(paulistring_int, bitindex)
    else
        paulistring_int = _setbittozero(paulistring_int, bitindex)
    end
    return paulistring_int
end

function _setbittoone(paulistring_int::Integer, bitindex::Integer)
    # set bit at bitindex to 1

    # shift ...00100...  to bitindex
    shifted_onebit = (typeof(paulistring_int)(1) << bitindex)

    # OR with operator to make sure that that bit is 1
    return paulistring_int | shifted_onebit
end

function _setbittozero(paulistring_int::Integer, bitindex::Integer)
    # set bit at bitindex to 0

    # flip all bits
    paulistring_int = ~paulistring_int

    # set target bit to one
    paulistring_int = _setbittoone(paulistring_int, bitindex)

    # flip all bits back, only the target bit is 0
    paulistring_int = ~paulistring_int
    return paulistring_int
end

function _readbit(oper::Integer, bitindex::Integer)
    # return integer with ...000[bit].

    # shift by bitindex
    shifted_oper = (oper >> bitindex)

    # AND with 1 to get first bit
    return shifted_oper & typeof(oper)(1)
end

@generated function alternatingmask(int::T) where {T<:Integer}
    # define our super bit mask looking like ....1010101.

    # length is the number of bits in the integer
    n_bits = min(bitsize(int), 2_048)  # for max 1024 qubits.
    mask = int(0)
    for ii in 0:(n_bits-1)
        if ii % 2 == 0
            mask = _setbittoone(mask, ii)
        end
    end
    return mask
end