# Defines mapping of integers 0, 1, 2, 3 to symbols :I, :X, :Y, :Z
const pauli_symbols::Vector{Symbol} = [:I, :X, :Y, :Z]

"""
    symboltoint(pstr::Union{Vector{Symbol}, Symbol})

Maps a symbol or a vector of symbols `pstr` to an integer Pauli string.

# Example
```
symboltoint([:X, :I])

# output

0x01
```
"""
function symboltoint(pstr)
    nqubits = length(pstr)
    converted_pstr = getinttype(nqubits)(0)
    for (qind, pauli) in enumerate(pstr)
        converted_pstr = setpauli(converted_pstr, pauli, qind)
    end
    return converted_pstr
end

"""
    symboltoint(nqubits::Integer, pstr::Vector{Symbol}, qinds::Vector{Int})

Maps a vector of symbols `pstr` acting on the indices `qinds` to an integer Pauli string. Other sites are set to the identity.
`qinds` can be any iterable.
"""
function symboltoint(nqubits::Integer, pstr, qinds)
    if nqubits < maximum(qinds)
        throw(ArgumentError("Indices in `qinds`=$qinds acts on more qubits than `nqubits`=$nqubits."))
    end
    inttype = getinttype(nqubits)
    return symboltoint(inttype, pstr, qinds)
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
    symboltoint(::TermType, pauli::Symbol, qind::Integer)

Maps a single symbol `pauli` acting on the index `qind` to an integer Pauli string with type `TermType`.
Other sites are set to the identity.
"""
function symboltoint(::Type{TT}, pauli::Symbol, qind::Integer) where {TT<:PauliStringType}
    converted_pauli = zero(TT)
    converted_pauli = setpauli(converted_pauli, pauli, qind)
    return converted_pauli
end

"""
    symboltoint(::TermType, pstr, qinds)

Maps a vector of symbols `pstr` acting on the indices `qinds` to an integer Pauli string with type `TermType`.
Other sites are set to the identity.
`qinds` can be any iterable.
"""
function symboltoint(::Type{TT}, pstr, qinds) where {TT<:PauliStringType}
    converted_pstr = zero(TT)
    for (qind, pauli) in zip(qinds, pstr)
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
symboltoint(pauli::Symbol) = findfirst(s -> s == pauli, pauli_symbols) - 1

"""
    inttosymbol(pauli::PauliType)

Maps an integer Pauli to its corresponding symbol.
"""
inttosymbol(pauli::PauliType) = pauli_symbols[pauli+1]

## get and set functions

"""
    getpauli(pstr::PauliStringType, index::Integer)

Gets the Pauli on index `index` of an integer Pauli string.
"""
function getpauli(pstr::PauliStringType, index::Integer)
    return _getpaulibits(pstr, index)
end


"""
    getpauli(pstr::PauliStringType, qinds::Vector{Integer})

Gets the Paulis on indices `qinds` of a `pstr` in the integer representation.
"""
function getpauli(pstr::PauliStringType, qinds)
    new_pstr = zero(pstr)
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
"""
inttostring(pstr::PauliType, nqubits::Integer) = prod("$(inttosymbol(getpauli(pstr, ii)))" for ii in 1:nqubits)

"""
Pretty string function.
"""
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
