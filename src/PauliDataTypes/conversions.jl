"""
    PauliSum(pstr::PauliString)

Constructor for a `PauliSum` on `nqubits` qubits from a `PauliString`.
"""
PauliSum(pstr::PauliString) = PauliSum(pstr.nqubits, pstr)

"""
    PauliSum(nq::Integer, pstr::PauliString)

Constructor for a `PauliSum` on `nqubits` qubits from a `PauliString`.
"""
function PauliSum(nq::Integer, pstr::PauliString{TT,CT}) where {TT,CT}
    _checknumberofqubits(nq, pstr)
    return PauliSum(nq, Dict{TT,CT}(pstr.term => pstr.coeff))
end

"""
    PauliSum(pstrs::Vector{PauliString})

Constructor for a `PauliSum` on `nqubits` qubits from a `PauliString`.
"""
function PauliSum(pstrs::Union{AbstractArray,Tuple,Base.Generator})
    nq = _checknumberofqubits(pstrs)

    CType = promote_type(coefftype.(pstrs)...)
    psum = PauliSum(CType, nq)

    sizehint!(psum, length(pstrs))

    for pstr in pstrs
        add!(psum, pstr)
    end
    return psum
end


PauliSum(vpsum::VectorPauliSum) = PauliSum(vpsum.nqubits, topaulistrings(vpsum))

VectorPauliSum(pstr::PauliString) = VectorPauliSum(pstr.nqubits, [pstr.term], [pstr.coeff])
VectorPauliSum(psum::PauliSum) = VectorPauliSum(psum.nqubits, collect(paulis(psum)), collect(coefficients(psum)))


"""
    tonumber(val::Number)

Trivial function returning a numerical value of a number.
Will be overloaded for custom wrapper types like `PathProperties`.
"""
function tonumber(val::Number)
    return val
end

# Conversion of PauliString and PauliSum to different coefficient types
function Base.convert(::Type{PauliString{TT1,CT1}}, pstr::PauliString{TT2,CT2}) where {TT1,TT2,CT1,CT2}
    if TT1 != TT2
        throw(ArgumentError("Cannot change term type from $TT2 to $TT1"))
    end
    return PauliString(pstr.nqubits, convert(TT1, pstr.term), convert(CT1, pstr.coeff))
end

function Base.convert(::Type{PauliSum{TT1,CT1}}, psum::PauliSum{TT2,CT2}) where {TT1,TT2,CT1,CT2}
    if TT1 != TT2
        throw(ArgumentError("Cannot change term type from $TT2 to $TT1"))
    end
    return PauliSum(psum.nqubits, convert(Dict{TT1,CT1}, psum.terms))

end

# # Examples
# ```julia
# convertcoefftype(Float64, PauliString(2, :X, 1, 1+0im))
# ```
# """
function convertcoefftype(::Type{CT1}, pstr::PauliString{TT,CT2}) where {TT,CT1,CT2}
    return PauliString(pstr.nqubits, pstr.term, convert(CT1, pstr.coeff))
end

# # Examples
# ```julia
# psum = PauliSum(PauliString(2, :X, 1, 1+0im))
# convertcoefftype(Float64, psum)
# ```
# """
function convertcoefftype(::Type{CT1}, psum::PauliSum{TT,CT2}) where {TT,CT1,CT2}
    return PauliSum(psum.nqubits, convert(Dict{TT,CT1}, psum.terms))
end



# Checks whether the number of qubits `nqubits` is the same between our datatypes.
function _checknumberofqubits(nq::Integer, pobj)
    obj_nq = nqubits(pobj)
    if nq != obj_nq
        throw(
            ArgumentError(
                "Number of qubits ($(nq)) must equal number of qubits ($(obj_nq)) in $(typeof(pobj))"
            )
        )
    end
    return nq
end



# Checks whether the number of qubits `nqubits` is the same between our datatypes.
function _checknumberofqubits(pobj1, pobj2)
    nq1 = nqubits(pobj1)
    nq2 = nqubits(pobj2)
    if nq1 != nq2
        throw(
            ArgumentError(
                "Number of qubits ($(nq1)) in $(typeof(pobj1)) must equal number of qubits ($(nq2)) in $(typeof(pobj2))"
            )
        )
    end
    return nq1
end

"""
Checks whether the number of qubits `nqubits` is the same between in some collection.
"""
function _checknumberofqubits(pobjects::Union{AbstractArray,Tuple,Base.Generator})

    if !allequal(nqubits(pobj) for pobj in pobjects)
        throw(
            ArgumentError(
                "Number of qubits in passed collection of type $(typeof(pobjects)) is not consistent."
            )
        )
    end
    return nqubits(first(pobjects))
end


# throw error for miss-matched term/Pauli types

function _checktermtype(pobj1, pobj2)
    if paulitype(pobj1) != paulitype(pobj2)
        throw(ArgumentError("Pauli types do not match. Got $(paulitype(pobj1)) and $(paulitype(pobj2))."))
    end
end