function countweight(pauli_op::Vector{Symbol}, args...; kwargs...)
    return sum(op in (:X, :Y, :Z) for op in pauli_op)
end

function countweight(oper::PauliStringType; kwargs...)
    return countbitweight(oper; kwargs...)
end


function countxy(oper::PauliStringType; kwargs...)
    return countbitxy(oper; kwargs...)
end

function countyz(oper::PauliStringType; kwargs...)
    return countbityz(oper; kwargs...)
end

containsXorY(symbs::AbstractArray{Symbol}, args...) = :X in symbs || :Y in symbs
containsXorY(int::PauliStringType, args...) = countxy(int) > 0
containsYorZ(int::PauliStringType, args...) = countyz(int) > 0


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

function commutes(sym1::Symbol, pauli_ind::SinglePauliType)::Bool
    return commutes(sym1, inttosymbol(pauli_ind))
end

function commutes(oper1::PauliStringType, oper2::PauliStringType)
    return bitcommutes(oper1, oper2)
end

function commutes(oper1::Dict{T,Float64}, oper2::Dict{T,Float64}) where {T<:PauliStringType}
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

function commutator(oper1::PauliStringType, oper2::PauliStringType)

    if commutes(oper1, oper2)
        total_sign = ComplexF64(0.0)
        new_oper = zero(typeof(oper1))
    else
        total_sign, new_oper = pauliprod(oper1, oper2)
    end
    # commutator is [A, B] = AB - BA = 2AB for non-commuting (meaning anti-commuting) Paulis
    return 2 * total_sign, new_oper
end

function commutator(op_dict1::Dict{OpType,CoeffType1}, op_dict2::Dict{OpType,CoeffType2}) where {OpType<:PauliStringType,CoeffType1,CoeffType2}
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
function pauliprod(pstr1::PauliString, pstr2::PauliString)
    _checknumberofqubits(pstr1, pstr2)
    sign, coeff = pauliprod(pstr1.operator, pstr2.operator)
    return PauliString(pstr1.nqubits, coeff, sign * pstr1.coeff * pstr2.coeff)
end

function pauliprod(op1::PauliStringType, op2::PauliStringType)
    # This function is for when we need to globally check the sign of the product (like in general products of Paulis, not local Pauli gates)
    op3 = bitpaulimultiply(op1, op2)
    sign = calculatesign(op1, op2, op3)
    return sign, op3
end

function pauliprod(op1::Symbol, op2::SinglePauliType)
    # assume that just one qubit is involved because we check commutation with a single Symbol
    return pauliprod(symboltoint(op1), op2, 1:1)
end

function pauliprod(op1::PauliStringType, op2::PauliStringType, changed_indices)
    # Calculate the Pauli product when you know on which sites the Paulis differ (changed_indices)
    op3 = bitpaulimultiply(op1, op2)
    sign = calculatesign(op1, op2, op3, changed_indices)
    return sign, op3
end

function calculatesign(op1::PauliStringType, op2::PauliStringType, op3::PauliStringType)
    # Calculate the sign of the product, loop as long as neither of the operators are Identity
    sign = Complex{Int64}(1)
    identity_pauli = 0
    while op1 > identity_pauli || op2 > identity_pauli  # while both are not identity
        sign *= calculatesign(op1, op2, op3, 1:1)
        op1 = _paulishiftright(op1)
        op2 = _paulishiftright(op2)
        op3 = _paulishiftright(op3)
    end
    return sign
end
function calculatesign(op1::PauliStringType, op2::PauliStringType, op3::PauliStringType, changed_indices)
    # Calculate the sign of the product but when you know on which sites the Paulis differ (changed_indices)
    sign = Complex{Int64}(1)
    for qind in changed_indices
        sign *= generalizedlevicivita(
            getpaulielement(op1, qind),
            getpaulielement(op2, qind),
            getpaulielement(op3, qind)
        )
    end
    return sign
end

const generalized_levicivita_matrix = permutedims(cat(
        [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1], # first arg is I
        [0 1 0 0; 1 0 0 0; 0 0 0 1im; 0 0 -1im 0], # first arg is X
        [0 0 1 0; 0 0 0 -1im; 1 0 0 0; 0 1im 0 0], # first arg is Y
        [0 0 0 1; 0 0 1im 0; 0 -1im 0 0; 1 0 0 0]; # first arg is Z
        dims=3), (2, 3, 1))

function generalizedlevicivita(n1::SinglePauliType, n2::SinglePauliType, n3::SinglePauliType)
    # acts like levicivita but yields the correct sign for products with I or P^2, and takes care of the imaginary coefficients in Pauli products
    return generalized_levicivita_matrix[n1+1, n2+1, n3+1]
end