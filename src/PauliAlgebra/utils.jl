
const pauli_symbols::Vector{Symbol} = [:I, :X, :Y, :Z]   # maps to 0 1 2 3

symboltoint(pauli::Symbol) = findfirst(s -> s == pauli, pauli_symbols) - 1
symboltoint(pauli) = pauli
inttosymbol(pauli::PauliType) = pauli_symbols[pauli+1]
inttosymbol(pauli) = pauli

function symboltoint(pstr::Vector{Symbol})
    nqubits = length(pstr)
    converted_pstr = getinttype(nqubits)(0)
    for (qind, pauli) in enumerate(pstr)
        converted_pstr = setpauli(converted_pstr, pauli, qind)
    end
    return converted_pstr
end

function symboltoint(nqubits::Integer, pstr::Vector{Symbol}, qinds) # TODO: pauli argument should always be first!
    inttype = getinttype(nqubits)
    converted_pstr = inttype(0)
    for (qind, pauli) in zip(qinds, pstr)
        converted_pstr = setpauli(converted_pstr, pauli, qind)
    end
    return converted_pstr
end

function symboltoint(nqubits::Integer, pauli::Symbol, qind::Integer)
    inttype = getinttype(nqubits)
    converted_pauli = inttype(0)
    converted_pauli = setpauli(converted_pauli, pauli, qind)
    return converted_pauli
end

function inttosymbol(pstr::PauliStringType, nqubits::Integer)
    converted_pstr = [:I for _ in 1:nqubits]
    for ii in 1:nqubits
        converted_pstr[ii] = inttosymbol(getpauli(pstr, ii))
    end
    return converted_pstr
end

## get and set functions
function getpauli(pstr::PauliString, index::Integer)
    return getpauli(pstr.operator, index)
end

function getpauli(pstr::PauliStringType, index::Integer)
    return _getpaulibits(pstr, index)
end

function setpauli(pstr::PauliString, pauli::T, index::Integer) where {T<:Union{Symbol,PauliType}}
    return PauliString(pstr.nqubits, setpauli(pstr.operator, pauli, index), str.coeff)
end

function setpauli(pstr::PauliStringType, pauli::PauliType, index::Integer)
    return _setpaulibits(pstr, pauli, index)
end

function setpauli(pstr::PauliStringType, pauli::Symbol, index::Integer)
    # `symboltoint` to ensure we work with `PauliType`, i.e., integers
    return setpauli(pstr, symboltoint(pauli), index)
end

## Helper functions for pretty printing
inttostring(pstr::PauliType, nqubits) = prod("$(inttosymbol(getpauli(pstr, ii)))" for ii in 1:nqubits)

function getprettystr(psum::Dict, nqubits::Int; max_lines=20)
    str = ""
    header = length(psum) == 1 ? "1 Pauli term: \n" : "$(length(psum)) Pauli terms:\n"
    str *= header

    for (ii, (op, coeff)) in enumerate(psum)
        if ii > max_lines
            new_str = "  ⋮"
            str *= new_str
            break
        end
        pauli_string = inttostring(op, nqubits)
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