# testing type aliases here
const PauliStringType = Union{UInt8,UInt16,UInt32,UInt64,UInt128,UInt256,BigInt,Int} # to be maintained when we adapt getinttype()
const SinglePauliType = PauliStringType


function getinttype(nq)
    nbits = 2 * nq
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
    elseif nbits <= 256    # These are super slow to hash
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

function countbitweight(oper::PauliStringType; kwargs...)
    mask = alternatingmask(oper)
    m1 = oper & mask
    m2 = oper & (mask << 1)
    res = m1 | (m2 >> 1)
    return count_ones(res)
end


function countbitxy(oper::PauliStringType; kwargs...)
    mask = alternatingmask(oper)

    op = oper ⊻ (oper >> 1)
    op = op & mask

    return count_ones(op)
end

function countbityz(oper::PauliStringType; kwargs...)
    mask = alternatingmask(oper)

    op = oper & (mask << 1)

    return count_ones(op)
end

function bitcommutes(op1::PauliStringType, op2::PauliStringType)

    mask0 = alternatingmask(op1)
    mask1 = mask0 << 1

    # obtain the left (then right) bits of each Pauli pair, in-place
    aBits0 = mask0 & op1
    aBits1 = mask1 & op1
    bBits0 = mask0 & op2
    bBits1 = mask1 & op2

    # shift left bits to align with right bits
    aBits1 = aBits1 >> 1
    bBits1 = bBits1 >> 1

    # sets '10' at every Pauli index where individual pairs don't commute
    flags = (aBits0 & bBits1) ⊻ (aBits1 & bBits0)

    # strings commute if parity of non-commuting pairs is even
    return (count_ones(flags) % 2) == 0
end

bitpaulimultiply(op1::PauliStringType, op2::PauliStringType) = op1 ⊻ op2

paulishiftright(op::PauliStringType) = op >> 2  # this will truncate the first encoded Pauli string

function getbitpaulielement(oper::PauliStringType, index::Integer)
    bitindex = 2 * (index - 1)
    # shift to the right and mask with 3 (00000011) to get the two bits
    shifted_oper = (oper >> bitindex)
    return shifted_oper & UInt8(3)
end


function setbitpaulielement!(oper::PauliStringType, index::Integer, element::SinglePauliType)
    bitindex = 2 * (index - 1)

    b1 = _readbit(element, 0)
    b2 = _readbit(element, 1)

    oper = _setbit(oper, bitindex, b1)
    oper = _setbit(oper, bitindex + 1, b2)
    return oper
end

function _setbit(oper::PauliStringType, bitindex::Integer, bit::Integer)
    if bit == true  # set to one
        oper = _setbittoone(oper, bitindex)
    else
        oper = _setbittozero(oper, bitindex)
    end
    return oper
end

function _setbittoone(oper::Integer, bitindex::Integer)
    return oper | (typeof(oper)(1) << bitindex)
end

function _setbittozero(oper::Integer, bitindex::Integer)
    oper = ~oper  # flip all bits
    oper = _setbittoone(oper, bitindex)
    oper = ~oper # flip all bits back
    return oper
end

function _readbit(oper::Integer, bitindex::Integer)
    # bitindex begins at 0
    return (oper >> bitindex) & UInt8(1)
end

@generated function alternatingmask(int::T) where {T<:Integer}
    n_bits = min(bitsize(int), 2_048)  # for 1024 qubits
    mask = int(0)
    for ii in 0:(n_bits-1)
        if ii % 2 == 0
            mask = _setbittoone(mask, ii)
        end
    end
    return mask
end