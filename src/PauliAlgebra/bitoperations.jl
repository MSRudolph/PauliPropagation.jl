"""
    getinttype(nqubits::Integer)

Function to return the smallest integer type that can hold `nqubits`.
This is the type that will be used internally for representing Pauli strings.
"""
function getinttype(nqubits::Integer)
    # we need 2 bits per qubit
    nbits = 2 * nqubits


    # just over 8.3 Million is the largest integer type we can generate
    for trial_bits in nbits:2:8_300_000

        # we can check if the number of bits is divisible by 8
        # othervise we know it cannot be defined
        if !(trial_bits % 8 == 0)
            continue
        end

        # special clauses for inbuilt integer types
        if trial_bits == 8
            return UInt8
        elseif trial_bits == 16
            return UInt16
        elseif trial_bits == 32
            return UInt32
        elseif trial_bits == 64
            return UInt64
        end
        # stop at 64 bits because I am suspicious of UInt128

        trial_inttype_expr = Symbol("UInt", trial_bits)
        # check if the integer type is defined to avoid overrides
        if isdefined(PauliPropagation, trial_inttype_expr)
            return eval(trial_inttype_expr)
        end

        # defining the integer type can fail for bit numbers that are odd not not natively supported
        # just try the next number if that happens
        try
            @eval @define_integers $trial_bits
            # return the newly defined unsigned integer type
            return eval(trial_inttype_expr)
        catch ErrorException
            continue
        end
    end

    # if we reach here, we have failed to define the integer type
    # Falling back to BigInt
    @warn "Failed to define integer types for $nqubits qubits. Falling back to BigInt."
    return BigInt
end


# A priated hash function for unsigned integers from BitIntegers.jl
# This hashes for the value of the integer and is a lot faster than the default hash function.
Base.hash(v::BitIntegers.AbstractBitUnsigned, h::UInt) = Base.hash_integer(v, h)


# This function counts the number of 00 bit pairs in the integer Pauli string.
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


# This function counts the number of 01 (X) or 10 (Y) bit pairs in the integer Pauli string.
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


# This function counts the number of 10 (Y) or 11 (Z) bit pairs in the integer Pauli string.
function _countbityz(pstr::PauliStringType)
    # we use that both have a 1 on the left bit

    # get our super bit mask looking like ....1010101.
    mask = alternatingmask(pstr)

    # AND with the shifted mask to extract the 1's on the left bit
    op = pstr & (mask << 1)

    # count 1's to get the number of Y or Z Paulis
    return count_ones(op)
end


# This function checks if two integer Pauli strings commute.
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


# XOR between two Pauli different non-identity strings gives the third one. Ignores signs or any coefficient.
_bitpaulimultiply(pstr1::PauliStringType, pstr2::PauliStringType) = pstr1 ⊻ pstr2



# Shift to the right and truncate the first encoded Pauli string. Just a utility function.
_paulishiftright(pstr::PauliStringType) = pstr >> 2


# Create a mask for n bits from index1 to index2: n = index2 - index1 + 1
# to get the first n bits or int repr (2^n -1) 
_maskfromindices(type, index1::Integer, index2::Integer) = 
    type(2^(2 * (index2 - index1 + 1)) - 1)

# This function extracts the Pauli at position `index` from the integer Pauli string.
function _getpaulibits(pstr::PauliStringType, index::Integer)
    # we need to shift the integer by 2 * (index - 1), then the first two bits are target Pauli
    bitindex = 2 * (index - 1)

    # shift to the right
    shifted_pstr = (pstr >> bitindex)

    # AND with 3 (00000011) to get the first two bits
    return shifted_pstr & typeof(pstr)(3)
end


# This function extracts the Pauli from `index1` to `index2`.
function _getpaulibits(pstr::PauliStringType, index1::Integer, index2::Integer)
    bitindex = 2 * (index1 - 1)

    # shift to the right
    shifted_pstr = (pstr >> bitindex)

    return shifted_pstr & _maskfromindices(typeof(pstr), index1, index2)
end


# Gets the bit at index `bitindex` in the integer Pauli string.
function _getbit(pauli::Integer, bitindex::Integer)
    # return integer with ...000[bit].

    # shift by bitindex
    shifted_pauli = (pauli >> bitindex)

    # AND with 1 to get first bit
    return shifted_pauli & typeof(pauli)(1)
end


# This function sets the Pauli at position `index` in the integer Pauli string to `target_pauli`.
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


# This function sets the Pauli from `index1` to `index2` to `target_pstr`.
function _setpaulibits(pstr::PauliStringType, target_pstr::PauliStringType, index1::Integer, index2::Integer)
   
    bitindex = 2 * (index1 - 1)

    shifted_mask =  (_maskfromindices(typeof(pstr), index1, index2) << bitindex)

    # set bits to target pstr
    pstr = (pstr & ~shifted_mask) | (target_pstr << bitindex)
    return pstr
end


# Sets a bit at index `bitindex` in the integer Pauli string to the value of `target_bit`.
function _setbit(pstr::PauliStringType, target_bit::Integer, bitindex::Integer)
    # set bit at bitindex to bit

    if target_bit == true  # set to one
        pstr = _setbittoone(pstr, bitindex)
    else
        pstr = _setbittozero(pstr, bitindex)
    end
    return pstr
end


# Sets a bit at index `bitindex` in the integer Pauli string to 1.
function _setbittoone(pstr::Integer, bitindex::Integer)
    # set bit at bitindex to 1

    # shift ...00100...  to bitindex
    shifted_onebit = (typeof(pstr)(1) << bitindex)

    # OR with pauli string to make sure that that bit is 1
    return pstr | shifted_onebit
end



# Sets a bit at index `bitindex` in the integer Pauli string to 0.
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


# This mask helps us to parallelize the bit operations over all qubits.
@generated function alternatingmask(pstr::T) where {T<:PauliStringType}
    # define our super bit mask looking like ....1010101.

    # length is the number of bits in the integer
    n_bits = min(bitsize(pstr), 2_048)  # for max 1024 qubits.
    mask = zero(pstr)
    for ii in 0:(n_bits-1)
        if ii % 2 == 0
            mask = _setbittoone(mask, ii)
        end
    end
    return mask
end