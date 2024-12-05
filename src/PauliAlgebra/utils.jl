# Defines mapping of integers 0, 1, 2, 3 to symbols :I, :X, :Y, :Z
const pauli_symbols::Vector{Symbol} = [:I, :X, :Y, :Z]

"""
    symboltoint(pstr)

Maps a vector of symbols `pstr` to an integer Pauli string.
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
    symboltoint(nqubits::Integer, pstr, qinds)

Maps a vector of symbols `pstr` acting on the indices `qinds` to an integer Pauli string. Other sites are set to the identity.
`qinds` can be any iterable.
"""
function symboltoint(nqubits::Integer, pstr, qinds)
    inttype = getinttype(nqubits)
    converted_pstr = inttype(0)
    for (qind, pauli) in zip(qinds, pstr)
        converted_pstr = setpauli(converted_pstr, pauli, qind)
    end
    return converted_pstr
end

"""
    symboltoint(nqubits::Integer, pauli::Symbol, qind::Integer)

Maps a single symbol `pauli` acting on the index `qind` to an integer Pauli string. Other sites are set to the identity.
"""
function symboltoint(nqubits::Integer, pauli::Symbol, qind::Integer)
    inttype = getinttype(nqubits)
    converted_pauli = inttype(0)
    converted_pauli = setpauli(converted_pauli, pauli, qind)
    return converted_pauli
end

"""
    inttosymbol(pstr::PauliStringType, nqubits::Integer)

Maps an integer Pauli string to a vector of symbols.
"""
function inttosymbol(pstr::PauliStringType, nqubits::Integer)  # TODO: does the argument order need to change?
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
    getpauli(pstr::PauliString, index::Integer)

Gets the Pauli on index `index` of a `PauliString` in the integer representation.
"""
function getpauli(pstr::PauliString, index::Integer)
    return getpauli(pstr.term, index)
end

"""
    getpauli(pstr::PauliStringType, index::Integer)

Gets the Pauli on index `index` of an integer Pauli string.
"""
function getpauli(pstr::PauliStringType, index::Integer)
    return _getpaulibits(pstr, index)
end

"""
    setpauli(pstr::PauliString, target_pauli::T, index::Integer) where {T<:Union{Symbol,PauliType}}

Sets the Pauli on index `index` of an integer `PauliString` to `target_pauli`. That Pauli can be provided as integer (0, 1, 2, 3) or as a symbol (:I, :X, :Y, :Z).
"""
function setpauli(pstr::PauliString, target_pauli::T, index::Integer) where {T<:Union{Symbol,PauliType}}
    return PauliString(pstr.nqubits, setpauli(pstr.term, target_pauli, index), str.coeff)
end

"""
    setpauli(pstr::PauliStringType, target_pauli::PauliType, index::Integer)

Sets the Pauli on index `index` of an integer Pauli string to `target_pauli`. That Pauli should be provided as integer (0, 1, 2, 3).
"""
function setpauli(pstr::PauliStringType, target_pauli::PauliType, index::Integer)
    return _setpaulibits(pstr, target_pauli, index)
end

"""
    setpauli(pstr::PauliStringType, target_pauli::PauliType, index::Integer)

Sets the Pauli on index `index` of an integer Pauli string to `target_pauli`. That Pauli should be provided as a symbol (:I, :X, :Y, :Z).
"""
function setpauli(pstr::PauliStringType, target_pauli::Symbol, index::Integer)
    # `symboltoint` to ensure we work with `PauliType`, i.e., integers
    return setpauli(pstr, symboltoint(target_pauli), index)
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
        elseif isa(coeff, PathProperties)
            if isa(coeff.coeff, Number)
                coeff_str = "PathProperty($(round(coeff.coeff, sigdigits=5)))"
            else
                coeff_str = "PathProperty($(typeof(coeff.coeff)))"
            end
        else
            coeff_str = "($(typeof(coeff)))"
        end
        new_str = " $(coeff_str) * $(pauli_string)\n"
        str *= new_str
    end

    return str

end