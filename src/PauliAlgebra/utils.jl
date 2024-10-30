
const pauli_ops::Vector{Symbol} = [:I, :X, :Y, :Z]   # maps to 0 1 2 3

symboltoint(sym::Symbol) = findfirst(s -> s == sym, pauli_ops) - 1
symboltoint(i::Integer) = i
inttosymbol(int::Integer) = pauli_ops[int+1]
inttosymbol(s::Symbol) = s

function symboltoint(oper::Vector{Symbol})
    nq = length(oper)
    intop = getinttype(nq)(0)
    for (ii, symb) in enumerate(oper)
        intop = setpaulielement!(intop, ii, symboltoint(symb))
    end
    return intop
end

function symboltoint(nq::Integer, symbols::Vector{Symbol}, qinds)
    inttype = getinttype(nq)
    intop = inttype(0)
    for (op, qind) in zip(symbols, qinds)
        intop = setpaulielement!(intop, qind, symboltoint(op))
    end
    return intop
end

function symboltoint(nq::Integer, symbol::Symbol, qind::Integer)
    inttype = getinttype(nq)
    intop = inttype(0)
    intop = setpaulielement!(intop, qind, symboltoint(symbol))
    return intop
end

function inttosymbol(int::Integer, n_qubits::Integer)
    symbs = [:I for _ in 1:n_qubits]
    for ii in 1:n_qubits
        symbs[ii] = inttosymbol(getpaulielement(int, ii))
    end
    return symbs
end

## Helper functions for pretty printing
inttostring(op::Integer, nq) = prod("$(inttosymbol(getpaulielement(op, ii)))" for ii in 1:nq)

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