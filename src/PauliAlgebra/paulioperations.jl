# TODO: generate these definitions with Macro's instead? Easier to maintain and less error-prone

function countweight(pstr::PauliString)
    return countweight(pstr.operator)
end

function countweight(psum::PauliSum)
    return countweight(psum.op_dict)
end

function countweight(pauli_dict::Dict{OpType,CoeffType}) where {OpType<:PauliStringType,CoeffType}
    return [countweight(pauli) for pauli in keys(pauli_dict)]
end

function countweight(pstr_int::PauliStringType)
    return _countbitweight(pstr_int)
end

function countxy(pstr::PauliString)
    return countxy(pstr.operator)
end

function countxy(psum::PauliSum)
    return countxy(psum.op_dict)
end

function countxy(pauli_dict::Dict{OpType,CoeffType}) where {OpType<:PauliStringType,CoeffType}
    return [countxy(pauli) for pauli in keys(pauli_dict)]
end

function countxy(pstr_int::PauliStringType)
    return _countbitxy(pstr_int)
end

function countyz(pstr::PauliString)
    return countyz(pstr.operator)
end

function countyz(psum::PauliSum)
    return countxy(psum.op_dict)
end

function countyz(pauli_dict::Dict{OpType,CoeffType}) where {OpType<:PauliStringType,CoeffType}
    return [countxy(pauli) for pauli in keys(pauli_dict)]
end

function countyz(oper::PauliStringType)
    return _countbityz(oper)
end

containsXorY(pstr::PauliString) = containsXorY(pstr.operator)
containsXorY(pstr_int::PauliStringType) = countxy(pstr_int) > 0
containsYorZ(pstr::PauliString) = containsYorZ(pstr.operator)
containsYorZ(pstr_int::PauliStringType) = countyz(pstr_int) > 0


### All the commutation check functions
function commutes(sym1::Symbol, sym2::Symbol)::Bool
    if sym1 == :I || sym2 == :I
        return true
    else
        return sym1 == sym2
    end
end

function commutes(sym1::Symbol, pauli_ind::PauliType)::Bool
    return commutes(sym1, inttosymbol(pauli_ind))
end

function commutes(oper1::PauliStringType, oper2::PauliStringType)
    return _bitcommutes(oper1, oper2)
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

function commutator(pstr1_int::PauliStringType, pstr2_int::PauliStringType)

    if commutes(pstr1_int, pstr2_int)
        total_sign = ComplexF64(0.0)
        new_oper = zero(typeof(pstr1_int))
    else
        total_sign, new_oper = pauliprod(pstr1_int, pstr2_int)
    end
    # commutator is [A, B] = AB - BA = 2AB for non-commuting (meaning anti-commuting) Paulis
    return 2 * total_sign, new_oper
end

function commutator(pauli_dict1::Dict{OpType,CoeffType1}, pauli_dict2::Dict{OpType,CoeffType2}) where {OpType<:PauliStringType,CoeffType1,CoeffType2}
    # different types of coefficients are allowed but not different types of operators

    new_pauli_dict = Dict{OpType,ComplexF64}()

    for (pauli1, coeff1) in pauli_dict1, (pauli2, coeff2) in pauli_dict2
        if !commutes(pauli1, pauli2)
            sign, new_op = commutator(pauli1, pauli2)
            new_pauli_dict[new_op] = get(new_pauli_dict, new_op, ComplexF64(0.0)) + sign * coeff1 * coeff2
        end
    end

    # Get rid of the pauli strings with zero coeffs
    # TODO: possibly combine this with the loop above
    for (k, v) in new_pauli_dict
        if abs(v) ≈ 0.0
            delete!(new_pauli_dict, k)
        end
    end

    return new_pauli_dict
end

## Pauli product
function pauliprod(pstr1::PauliString, pstr2::PauliString)
    _checknumberofqubits(pstr1, pstr2)
    sign, coeff = pauliprod(pstr1.operator, pstr2.operator)
    return PauliString(pstr1.nqubits, coeff, sign * pstr1.coeff * pstr2.coeff)
end

function pauliprod(pauli1::PauliStringType, pauli2::PauliStringType)
    # This function is for when we need to globally check the sign of the product (like in general products of Paulis, not local Pauli gates)
    pauli3 = _bitpaulimultiply(pauli1, pauli2)
    sign = calculatesign(pauli1, pauli2, pauli3)
    return sign, pauli3
end

function pauliprod(pauli1::Symbol, pauli2::PauliType)
    # assume that just one qubit is involved because we check commutation with a single Symbol
    return pauliprod(symboltoint(pauli1), pauli2, 1:1)
end

function pauliprod(pauli1::PauliStringType, pauli2::PauliStringType, changed_indices)
    # Calculate the Pauli product when you know on which sites the Paulis differ (changed_indices)
    pauli3 = _bitpaulimultiply(pauli1, pauli2)
    sign = calculatesign(pauli1, pauli2, pauli3, changed_indices)
    return sign, pauli3
end

function calculatesign(pauli1::PauliStringType, pauli2::PauliStringType, pauli3::PauliStringType)
    # Calculate the sign of the product, loop as long as neither of the operators are Identity
    sign = Complex{Int64}(1)
    identity_pauli = 0
    while pauli1 > identity_pauli || pauli2 > identity_pauli  # while both are not identity
        sign *= calculatesign(pauli1, pauli2, pauli3, 1:1)
        pauli1 = _paulishiftright(pauli1)
        pauli2 = _paulishiftright(pauli2)
        pauli3 = _paulishiftright(pauli3)
    end
    return sign
end
function calculatesign(pauli1::PauliStringType, pauli2::PauliStringType, pauli3::PauliStringType, changed_indices)
    # Calculate the sign of the product but when you know on which sites the Paulis differ (changed_indices)
    sign = Complex{Int64}(1)
    for qind in changed_indices
        sign *= generalizedlevicivita(
            getpauli(pauli1, qind),
            getpauli(pauli2, qind),
            getpauli(pauli3, qind)
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

function generalizedlevicivita(pauli1::PauliType, pauli2::PauliType, pauli3::PauliType)
    # acts like levicivita but yields the correct sign for products with I or P^2, and takes care of the imaginary coefficients in Pauli products
    return generalized_levicivita_matrix[pauli1+1, pauli2+1, pauli3+1]
end