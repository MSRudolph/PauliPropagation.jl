
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

function commutes(oper1::Integer, oper2::Integer)
    return bitcommutes(oper1, oper2)
end

## Commutator
# TODO: add here

### Code for the product of pauli_ops
function pauliprod(pstr1::PauliString, pstr2::PauliString)
    sign, new_operator = pauliprod(pstr1.operator, pstr2.operator, pstr1.nqubits)
    return PauliString(pstr1.nqubits, new_operator, pstr1.coeff * pstr2.coeff * sign)
end

function pauliprod(op1::Integer, op2::Integer, nq::Int)
    return pauliprod(op1, op2, 1:nq)
end

function pauliprod(op1::Symbol, op2::Integer) # assume that just one qubit is involved, TODO: This must change to one bit operation on all qubits
    return pauliprod(symboltoint(op1), op2, 1:1)
end

function singlepauliprod(op1::Integer, op2::Integer)
    return pauliprod(getelement(op1, 1), getelement(op2, 1), 1:1)
end

function pauliprod(op1::Integer, op2::Integer, changed_indices) # TODO: find a way to circumvent changed_indices being passed
    op3 = bitpauliprod(op1, op2)
    sign = calculatesign(op1, op2, op3, changed_indices)
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

const generalized_levicivita_matrix = cat(
    [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1], # first arg is I
    [0 1im 0 0; 1 0 0 0; 0 0 0 1im; 0 0 -1im 0], # first arg is X
    [0 0 1im 0; 0 0 0 -1im; 1 0 0 0; 0 1im 0 0], # first arg is Y
    [0 0 0 1im; 0 0 1im 0; 0 -1im 0 0; 1 0 0 0]; # first arg is Z
    dims=3)

function generalizedlevicivita(n1::Integer, n2::Integer, n3::Integer)
    # acts like levicivita but yields the correct sign for products with I or P^2
    return generalized_levicivita_matrix[n1+1, n2+1, n3+1]
end