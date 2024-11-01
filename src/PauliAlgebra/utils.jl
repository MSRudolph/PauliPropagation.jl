
const pauli_ops::Vector{Symbol} = [:I, :X, :Y, :Z]   # maps to 0 1 2 3

symboltoint(pauli::Symbol) = findfirst(s -> s == pauli, pauli_ops) - 1
symboltoint(pauli) = pauli
inttosymbol(pauli::PauliType) = pauli_ops[pauli+1]
inttosymbol(pauli) = pauli

function symboltoint(pstr_vec::Vector{Symbol})
    nqubits = length(pstr_vec)
    pstr_int = getinttype(nqubits)(0)
    for (ii, symbol) in enumerate(pstr_vec)
        pstr_int = setpauli(pstr_int, ii, symboltoint(symbol))
    end
    return pstr_int
end

function symboltoint(nqubits::Integer, pstr_vec::Vector{Symbol}, qinds) # TODO: pauli argument should always be first!
    inttype = getinttype(nqubits)
    pstr_int = inttype(0)
    for (pauli, qind) in zip(pstr_vec, qinds)
        pstr_int = setpauli(pstr_int, qind, symboltoint(pauli))
    end
    return pstr_int
end

function symboltoint(nqubits::Integer, pauli::Symbol, qind::Integer)
    inttype = getinttype(nqubits)
    pstr_int = inttype(0)
    pstr_int = setpauli(pstr_int, qind, symboltoint(pauli))
    return pstr_int
end

function inttosymbol(pauli::PauliStringType, nqubits::Integer)
    symbols = [:I for _ in 1:nqubits]
    for ii in 1:nqubits
        symbols[ii] = inttosymbol(getpauli(pauli, ii))
    end
    return symbols
end

## get and set functions
function getpauli(pstr::PauliString, index::Integer)
    return getpauli(pstr.operator, index)
end

function getpauli(pstr_int::PauliStringType, index::Integer)
    return getpaulibits(pstr_int, index)
end

function setpauli(pstr::PauliString, index::Integer, pauli::T) where {T<:Union{Symbol,PauliType}}
    return PauliString(pstr.nqubits, setpauli(pstr.operator, index, pauli), str.coeff)
end

function setpauli(pstr_int::PauliStringType, index, pauli::PauliType)
    return setpaulibits(pstr_int, index, pauli)
end

function setpauli(pstr_int::PauliStringType, index, pauli::Symbol)
    return setpauli(pstr_int, index, symboltoint(pauli))
end


## Helper functions for pretty printing
inttostring(pstr_int::PauliType, nq) = prod("$(inttosymbol(getpauli(pstr_int, ii)))" for ii in 1:nq)

function getprettystr(d::Dict, nq::Int; max_lines=20)
    str = ""
    header = length(d) == 1 ? "1 Pauli term: \n" : "$(length(d)) Pauli terms:\n"
    str *= header

    for (ii, (op, coeff)) in enumerate(d)
        if ii > max_lines
            new_str = "  ⋮"
            str *= new_str
            break
        end
        pauli_string = inttostring(op, nq)
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