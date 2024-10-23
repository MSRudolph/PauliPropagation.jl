
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

function commutes(oper1::Dict{T,Float64}, oper2::Dict{T,Float64}) where {T<:Integer}
    comm = commutator(oper1, oper2)
    return isempty(comm)
end

## Commutator
function commutator(psum1::PauliSum, psum2::PauliSum)
    new_op_dict = commutator(psum1.op_dict, psum2.op_dict)
    return PauliSum(psum1.nqubits, new_op_dict)
end

function commutator(pstr1::PauliString, pstr2::PauliString)
    new_coeff, new_op = commutator(pstr1.operator, pstr2.operator)
    return PauliString(pstr1.nqubits, new_op, new_coeff)
end

commutator(psum::PauliSum, pstr::PauliString) = commutator(psum, PauliSum(pstr))
commutator(pstr::PauliString, psum::PauliSum) = commutator(PauliSum(pstr), psum)

function commutator(oper1::T, oper2::T) where {T<:Integer}
    new_oper = zero(T)

    if commutes(oper1, oper2)
        total_sign = ComplexF64(0.0)
    else
        total_sign, new_oper = pauliprod(oper1, oper2)
    end
    return 2 * total_sign, new_oper
end

function commutator(op_dict1::Dict{OpType,CoeffType1}, op_dict2::Dict{OpType,CoeffType2}) where {OpType,CoeffType1,CoeffType2}
    # different types of coefficients are allowed but not different types of operators

    new_op_dict = Dict{OpType,ComplexF64}()

    for (op1, coeff1) in op_dict1, (op2, coeff2) in op_dict2
        if !commutes(op1, op2)
            sign, new_op = commutator(op1, op2)
            new_op_dict[new_op] = get(new_op_dict, new_op, ComplexF64(0.0)) + sign * coeff1 * coeff2
        end
    end

    # Get rid of the pauli strings with zero coeffs
    for (k, v) in new_op_dict
        if abs(v) ≈ 0.0
            delete!(new_op_dict, k)
        end
    end

    return new_op_dict
end

## Pauli product
function pauliprod(op1::Integer, op2::Integer)
    # TODO: Can this be done with a few bit operations?
    sign = Complex{Int64}(1)
    new_op = bitpauliprod(op1, op2)
    op3 = new_op

    # Calculate the sign of the product, loop as long as neither of the operators are Identity
    while op1 > 0 || op2 > 0
        sign *= PauliPropagation.calculatesign(op1, op2, op3, 1:1)
        op1 = bitshiftright(op1)
        op2 = bitshiftright(op2)
        op3 = bitshiftright(op3)
    end
    return sign, new_op
end

function pauliprod(op1::Symbol, op2::Integer) # assume that just one qubit is involved
    return pauliprod(symboltoint(op1), op2, 1:1)
end

function pauliprod(op1::Integer, op2::Integer, changed_indices)
    op3 = bitpauliprod(op1, op2)
    sign = calculatesign(op1, op2, op3, changed_indices)
    return sign, op3
end

function calculatesign(op1::Integer, op2::Integer, op3::Integer, changed_indices)
    sign = Complex{Int64}(1)
    for qind in changed_indices
        sign *= generalizedlevicivita(
            getelement(op1, qind),
            getelement(op2, qind),
            getelement(op3, qind)
        )
    end
    return sign
end

const generalized_levicivita_matrix = permutedims(cat(
    [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1], # first arg is I
    [0 1 0 0; 1 0 0 0; 0 0 0 1im; 0 0 -1im 0], # first arg is X
    [0 0 1 0; 0 0 0 -1im; 1 0 0 0; 0 1im 0 0], # first arg is Y
    [0 0 0 1; 0 0 1im 0; 0 -1im 0 0; 1 0 0 0]; # first arg is Z
    dims=3), (2,3,1))

function generalizedlevicivita(n1::Integer, n2::Integer, n3::Integer)
    # acts like levicivita but yields the correct sign for products with I or P^2
    return generalized_levicivita_matrix[n1+1, n2+1, n3+1]
end