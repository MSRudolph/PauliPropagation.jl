
const pauli_ops::Vector{Symbol} = [:I, :X, :Y, :Z]   # maps to 0 1 2 3

symboltoint(sym::Symbol) = findfirst(s -> s == sym, pauli_ops) - 1
symboltoint(i::Integer) = i
inttosymbol(int::Integer) = pauli_ops[int+1]
inttosymbol(s::Symbol) = s

function symboltoint(oper::Vector{Symbol})
    nq = length(oper)
    intoper = getinttype(nq)(0)
    for (ii, symb) in enumerate(oper)
        intoper = setelement!(intoper, ii, symboltoint(symb))
    end
    return intoper
end

function inttosymbol(int::Integer, n_qubits::Integer)
    symbs = [:I for _ in 1:n_qubits]
    for ii in 1:n_qubits
        symbs[ii] = inttosymbol(getelement(int, ii))
    end
    return symbs
end


inttostring(op::Unsigned) = prod("$(inttosymbol(getelement(op, ii)))" for ii in 1:Int(bitsize(op) / 2))
inttostring(op::Unsigned, nq) = prod("$(inttosymbol(getelement(op, ii)))" for ii in 1:nq)

import Base.show
function show(op::Integer)
    println(inttostring(op))

end

function show(op::Integer, n::Int)
    max_qubits_in_integer = round(Int, bitsize(typeof(op)) / 2)
    nind = min(max_qubits_in_integer, n)

    print_string = inttostring(op)[1:nind]
    println(print_string)

end


function show(d::Dict; max_lines=20)
    show(d, Int(bitsize(first(d)[1]) / 2); max_lines=max_lines)

end

function show(d::Dict, nq::Int; max_lines=20)
    println(getdictstr(d, nq; max_lines=max_lines))

end

function getdictstr(d::Dict, nq::Int; max_lines=20)
    str = ""
    header = "$(typeof(d)) with " * (length(d) == 1 ? "1 entry:\n" : "$(length(d)) entries:\n")
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
        new_str = "  $(pauli_string) => $coeff \n"
        str *= new_str
    end

    return str

end

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