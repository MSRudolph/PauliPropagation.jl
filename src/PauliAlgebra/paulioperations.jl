
function countweight(pauli_op::Vector{Symbol}, args...; kwargs...)
    return sum(op in (:X, :Y, :Z) for op in pauli_op)
end

function countweight(oper::Integer; kwargs...)
    return countbitweight(oper; kwargs...)
end


function countxy(oper::Integer; kwargs...)
    return countbitxy(oper; kwargs...)
end

function countyz(oper::Integer; kwargs...)
    return countbityz(oper; kwargs...)
end

containsXorY(symbs::AbstractArray{Symbol}, args...) = :X in symbs || :Y in symbs
containsXorY(int::Integer, args...) = countxy(int) > 0
containsYorZ(int::Integer, args...) = countyz(int) > 0


### All the commutation check functions

# TODO: Should we even support operations on non-integer operators?
function commutes(gate_generator::AbstractArray{T}, pauli_op::AbstractArray{T})::Bool where {T}
    return sum(!commutes(o1, o2) for (o1, o2) in zip(gate_generator, pauli_op)) % 2 == 0
end

function commutes(sym1::Symbol, sym2::Symbol)::Bool
    if sym1 == :I || sym2 == :I
        return true
    else
        return sym1 == sym2
    end
end

function commutes(sym1::Symbol, pauli_ind::Integer)::Bool
    return commutes(sym1, inttosymbol(pauli_ind))
end

function commutes(gate::PauliGateUnion, oper)
    return sum(!commutes(gate_sym, getelement(oper, qind)) for (qind, gate_sym) in zip(gate.qinds, gate.symbols)) % 2 == 0
end

function commutes(gate::FastPauliGate, oper::Integer)
    return commutes(gate.bitoperator, oper)
end

function commutes(oper1::Integer, oper2::Integer)
    return bitcommutes(oper1, oper2)
end

### Code for the product of pauli_ops


function pauliprod(pstr1::PauliString, pstr2::PauliString)
    sign, new_operator = pauliprod(pstr1.operator, pstr2.operator, pstr1.nqubits)
    return PauliString(pstr1.nqubits, new_operator, pstr1.coeff * pstr2.coeff * sign)
end

function pauliprod(op1::Integer, op2::Integer, nq::Int)
    return pauliprod(op1, op2, 1:nq)
end

function pauliprod(op1::Integer, op2::Integer, changed_indices) # TODO: find a way to circumvent changed_indices being passed
    op3 = bitpauliprod(op1, op2)
    sign = calculatesign(op1, op2, op3, changed_indices)
    return sign, op3
end


function pauliprod(op1::Symbol, op2::Integer) # assume that just one qubit is involved, TODO: This must change to one bit operation on all qubits
    return pauliprod(symboltoint(op1), op2, 1:1)
end

function singlepauliprod(op1::Integer, op2::Integer)
    op3 = getelement(op1, 1) ⊻ getelement(op2, 1)
    sign = calculatesign(op1, op2, op3, 1:1)
    return sign, op3
end

function calculatesign(op1::Integer, op2::Integer, op3::Integer, changed_indices)
    sign = 1
    for qind in changed_indices
        sign *= generalizedlevicivita(
            getelement(op1, qind),
            getelement(op2, qind),
            getelement(op3, qind)
        )
    end
    return sign
end


# function paulimult(sym1::Symbol, sym2::Symbol)
#     ind1 = symboltoint(sym1)
#     ind2 = symboltoint(sym2)
#     sign, ind3 = paulimult(ind1, ind2)
#     return sign, inttosymbol(ind3)
# end

# function paulimult(sym::Symbol, pauli_ind::Integer)
#     ind = symboltoint(sym)
#     sign, ind3 = paulimult(ind, pauli_ind)
#     return sign, ind3
# end

# function paulimult(pauli1::Integer, pauli2::Integer)
#     if pauli1 == 0
#         return 1, pauli2
#     elseif pauli2 == 0
#         return 1, pauli1
#     elseif pauli1 == pauli2
#         return 1, 0
#     else
#         new_pauli = 0
#         for ii in 1:3
#             if ii != pauli1 && ii != pauli2
#                 new_pauli = ii
#                 break
#             end
#         end
#         sign = levicivita(pauli1, pauli2, new_pauli)
#         return sign, new_pauli
#     end
# end

# function paulimult(op1::AbstractArray{T}, op2::AbstractArray{T}) where {T}  # TODO: should we even support this?
#     total_sign = -1
#     new_op = [:I for _ in eachindex(op1)]
#     for (ii, (o1, o2)) in enumerate(zip(op1, op2))
#         sign, new_o = paulimult(o1, o2)
#         total_sign *= sign
#         new_op[ii] = new_o
#     end
#     return total_sign, new_op
# end

# const levicivita_matrix = cat([0 0 0; 0 0 1; 0 -1 0],
#     [0 0 -1; 0 0 0; 1 0 0],
#     [0 1 0; -1 0 0; 0 0 0];
#     dims=3)

# function levicivita(n1::Integer, n2::Integer, n3::Integer)
#     return levicivita_matrix[n1, n2, n3]
# end


const generalized_levicivita_matrix = cat(
    [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1], # first arg is I
    [0 1 0 0; 1 0 0 0; 0 0 0 1; 0 0 -1 0], # first arg is X
    [0 0 1 0; 0 0 0 -1; 1 0 0 0; 0 1 0 0], # first arg is Y
    [0 0 0 1; 0 0 1 0; 0 -1 0 0; 1 0 0 0]; # first arg is Z
    dims=3)

function generalizedlevicivita(n1::Integer, n2::Integer, n3::Integer)
    # acts like levicivita but yields the correct sign for products with I or P^2
    return generalized_levicivita_matrix[n1+1, n2+1, n3+1]
end