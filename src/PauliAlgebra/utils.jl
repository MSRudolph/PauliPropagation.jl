# Defines mapping of integers 0, 1, 2, 3 to symbols :I, :X, :Y, :Z
const pauli_symbols = (:I, :X, :Y, :Z)


# Returns the integer representation of the identity Pauli string acting on `nqubits` qubits.
# The type of will be the smallest integer type that can hold the number of qubits, as given by `getinttype(nqubits)`.
function identitypauli(nqubits::Integer)
    return identitypauli(getinttype(nqubits))
end


# Returns the integer representation of the identity Pauli string with type `PauliStringType`.
function identitypauli(::Type{TT}) where {TT<:PauliStringType}
    return zero(TT)
end


# Returns an integer Pauli string of the same type as `pstr` with all Paulis set to identity.
function identitylike(::TT) where {TT<:PauliStringType}
    return identitypauli(TT)
end


"""
    symboltoint(pstr::Union{Vector{Symbol}, Symbol})

Maps a symbol or a vector of symbols `pstr` to an integer Pauli string.

Example:
```
symboltoint([:X, :I])
>>> 0x01
```
"""
function symboltoint(pstr)
    nqubits = length(pstr)
    converted_pstr = identitypauli(nqubits)
    for (qind, pauli) in enumerate(pstr)
        converted_pstr = setpauli(converted_pstr, pauli, qind)
    end
    return converted_pstr
end

"""
    symboltoint(nqubits::Integer, paulis::Vector{Symbol}, qinds::Vector{Int})

Maps a vector of symbols `pstr` acting on the indices `qinds` to an integer Pauli string. Other sites are set to the identity.
`qinds` can be any iterable.
"""
function symboltoint(nqubits::Integer, paulis, qinds)
    # whatever type these are, they should be the same length
    if length(paulis) != length(qinds)
        throw(ArgumentError("Length of `paulis=$(length(paulis))` and `qinds`=$(length(qinds)) should be the same."))
    end

    # check that qinds are in the correct range 1 <= qind <= nqubits
    if any(qind -> !(1 <= qind <= nqubits), qinds)
        throw(ArgumentError("Indices `qinds` should be in the range 1 <= ... <= nqubits=$nqubits. Got `qinds`=$(qinds)."))
    end

    # check that indices are unique, which is otherwise likely unintended
    if length(qinds) != length(Set(qinds))
        throw(ArgumentError("Indices `qinds` should be unique. Got `qinds`=$(qinds)."))
    end

    TT = getinttype(nqubits)
    return symboltoint(TT, paulis, qinds)
end

"""
    symboltoint(nqubits::Integer, pauli::Symbol, qind::Integer)

Maps a single symbol `pauli` acting on the index `qind` to an integer Pauli string. Other sites are set to the identity.
"""
function symboltoint(nqubits::Integer, pauli::Symbol, qind::Integer)
    TT = getinttype(nqubits)
    return symboltoint(TT, pauli, qind)
end

"""
    symboltoint(::PauliStringType, pauli::Symbol, qind::Integer)

Maps a single symbol `pauli` acting on the index `qind` to an integer Pauli string with type `PauliStringType`.
Other sites are set to the identity.
"""
function symboltoint(::Type{TT}, pauli::Symbol, qind::Integer) where {TT<:PauliStringType}
    converted_pauli = identitypauli(TT)
    converted_pauli = setpauli(converted_pauli, pauli, qind)
    return converted_pauli
end

"""
    symboltoint(::PauliStringType, paulis, qinds)

Maps a vector of symbols `paulis` acting on the indices `qinds` to an integer Pauli string with type `PauliStringType`.
Other sites are set to the identity.
`qinds` can be any iterable.
"""
function symboltoint(::Type{TT}, paulis, qinds) where {TT<:PauliStringType}
    converted_pstr = identitypauli(TT)
    for (qind, pauli) in zip(qinds, paulis)
        converted_pstr = setpauli(converted_pstr, pauli, qind)
    end
    return converted_pstr
end

"""
    inttosymbol(pstr::PauliStringType, nqubits::Integer)

Maps an integer Pauli string to a vector of symbols.
"""
function inttosymbol(pstr::PauliStringType, nqubits::Integer)
    converted_pstr = [:I for _ in 1:nqubits]
    for ii in 1:nqubits
        converted_pstr[ii] = inttosymbol(getpauli(pstr, ii))
    end
    return converted_pstr
end

"""
    symboltoint(pauli::Symbol)

Maps a single symbol to its corresponding integer representation.
"""
function symboltoint(pauli::Symbol)

    ind = findfirst(s -> s == pauli, pauli_symbols)
    if isnothing(ind)
        throw(ArgumentError("Symbol $pauli is not a valid Pauli symbol."))
    end
    return ind - 1
end

"""
    inttosymbol(pauli::PauliType)

Maps an integer Pauli to its corresponding symbol.
"""
function inttosymbol(pauli::PauliType)
    if !(0 <= pauli <= 3)
        throw(ArgumentError("Pauli $pauli is not a valid Pauli integer."))
    end
    return pauli_symbols[pauli+1]
end

# trivial functions that return the if it is already in the correct type
symboltoint(pauli::PauliStringType) = pauli
# TODO: we might want versions for tuples, arrays and vectors
inttosymbol(pauli::Symbol) = pauli


## testing for equality between integer and symbol representations
"""
    ispauli(pauli1::Union{Symbol, PauliType}, pauli2::Union{Symbol, PauliType})
    
    ispauli(pauli1::Union{Vector{Symbol}, PauliStringType}, pauli2::Union{Vector{Symbol}, PauliStringType})

Check if two Paulis are equal, where one is given as a symbol and the other as an integer.
"""
function ispauli(pauli1, pauli2)
    # always convert to integer representation
    return symboltoint(pauli1) == symboltoint(pauli2)
end


## get and set functions


# Gets the Pauli on index `index` of an integer Pauli string.
function getpauli(pstr::PauliStringType, index::Integer)
    return _getpaulibits(pstr, index)
end


"""
    getpauli(pstr::PauliStringType, qinds::Vector{Integer})

Gets the Paulis on indices `qinds` of a `pstr` in the integer representation.
"""
function getpauli(pstr::PauliStringType, qinds)
    new_pstr = identitylike(pstr)
    # Get the Paulis on the indices `qinds`
    for (ii, qind) in enumerate(qinds)
        pauli = getpauli(pstr, qind)
        new_pstr = setpauli(new_pstr, pauli, ii)
    end
    return new_pstr

end

"""
    setpauli(pstr::PauliStringType, target_pauli::PauliType, index::Integer)

Sets the Pauli on index `index` of an integer Pauli string to `target_pauli`. 
That Pauli should be provided as integer (0, 1, 2, 3).
"""
function setpauli(pstr::PauliStringType, target_pauli::PauliType, index::Integer)
    return _setpaulibits(pstr, target_pauli, index)
end

"""
    setpauli(pstr::PauliStringType, target_pauli::Symbol, index::Integer)

Sets the Pauli on `index` of an integer Pauli string to `target_pauli`. 
That Pauli should be provided as a symbol (:I, :X, :Y, :Z).
"""
function setpauli(pstr::PauliStringType, target_pauli::Symbol, index::Integer)
    # `symboltoint` to ensure we work with `PauliType`, i.e., integers
    return setpauli(pstr, symboltoint(target_pauli), index)
end


"""
    setpauli(
        pstr::PauliStringType, 
        target_paulis::PauliStringType, 
        qinds::Vector{Integer}
    )

Set the Paulis `qinds` of an integer Pauli string `pstr` to `target_paulis`.
Use Tuples for `qinds` in performance critical functions because they are immutable.
"""
function setpauli(pstr::PauliStringType, target_paulis::PauliStringType, qinds)
    for (ii, qind) in enumerate(qinds)
        pstr = setpauli(pstr, getpauli(target_paulis, ii), qind)
    end
    return pstr
end

"""
    setpauli(
        pstr::PauliStringType, 
        target_paulis::Vector{Symbol}, 
        qinds::Vector{Integer}
    )

Set the Paulis `qinds` of an integer Pauli string `pstr` to `target_paulis`.
`target_paulis` is a vector of symbols.
Use tuples in performance critical functions because they are immutable.
"""
function setpauli(pstr::PauliStringType, target_paulis, qinds)
    for (ii, qind) in enumerate(qinds)
        pstr = setpauli(pstr, symboltoint(target_paulis[ii]), qind)
    end
    return pstr
end

## Helper functions for pretty printing
"""
    inttostring(pstr::PauliType, nqubits::Integer)

Returns a string representation of an integer Pauli string `pstr` on `nqubits` qubits.
The characters of the string from left to right are the Paulis on the qubits from 1 to `nqubits`.
"""
inttostring(pstr::PauliType, nqubits::Integer) = prod("$(inttosymbol(getpauli(pstr, ii)))" for ii in 1:nqubits)


# Pretty string function.
function _getprettystr(psum::Dict, nqubits::Int; max_lines=20)
    # TODO: rework this pretty print to not build the string but keep streaming via show(io, ...)
    str = ""
    header = length(psum) == 1 ? "1 Pauli term: \n" : "$(length(psum)) Pauli terms:\n"
    str *= header

    for (ii, (pstr, coeff)) in enumerate(psum)
        if ii > max_lines
            new_str = "  ⋮"
            str *= new_str
            break
        end
        pauli_string = inttostring(pstr, nqubits)
        if length(pauli_string) > 20
            pauli_string = pauli_string[1:20] * "..."
        end
        if isa(coeff, Number)
            coeff_str = round(coeff, sigdigits=5)
        elseif isa(coeff, PathProperties) && hasfield(typeof(coeff), :coeff)
            PProp = string(typeof(coeff).name.name)
            if isa(coeff.coeff, Number)
                coeff_str = "$PProp($(round(coeff.coeff, sigdigits=5)))"
            else
                coeff_str = "$PProp($(typeof(coeff.coeff)))"
            end
        else
            coeff_str = "$(typeof(coeff))"
        end
        new_str = " $(coeff_str) * $(pauli_string)\n"
        str *= new_str
    end

    return str

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
    return PauliSum(pstr.nqubits, convert(Dict{TT1,CT1}, psum.terms))

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
