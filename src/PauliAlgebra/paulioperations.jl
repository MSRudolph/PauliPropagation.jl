# TODO: generate these definitions with Macro's instead? Easier to maintain and less error-prone

"""
    countweight(pstr::PauliString)

Function to count the weight of a `PauliString`.
"""
function countweight(pstr::PauliString)
    return countweight(pstr.operator)
end

"""
    countweight(pstr::PauliStringType)

Function to count the weight of an integer Pauli string.
"""
function countweight(pstr::PauliStringType)
    return _countbitweight(pstr)
end

"""
    countweight(psum::PauliSum)

Function to count the weight Pauli strings in a `PauliSum`. Returns an array of weights.
"""
function countweight(psum::PauliSum)
    return countweight(psum.op_dict)
end

"""
    countweight(psum::Dict)

Function to count the weight Pauli strings in a Pauli sum in the form of a `Dict`. Returns an array of weights.
"""
function countweight(psum::Dict)
    return [countweight(pstr) for pstr in keys(psum)]
end

"""
    countxy(pstr::PauliString)

Function to count the number of X and Y Paulis in a `PauliString`.
"""
function countxy(pstr::PauliString)
    return countxy(pstr.operator)
end

"""
    countxy(pstr::PauliStringType)

Function to count the number of X and Y Paulis in an integer Pauli string.
"""
function countxy(pstr::PauliStringType)
    return _countbitxy(pstr)
end

"""
    countxy(psum::PauliSum)

Function to count the number of X and Y Paulis in a `PauliSum`. Returns an array of counts.
"""
function countxy(psum::PauliSum)
    return countxy(psum.op_dict)
end

"""
    countxy(psum::Dict)

Function to count the number of X and Y Paulis in a Pauli sum in the form of a `Dict`. Returns an array of counts.
"""
function countxy(psum::Dict)
    return [countxy(pstr) for pstr in keys(psum)]
end

"""
    countyz(pstr::PauliString)

Function to count the number of Y and Z Paulis in a `PauliString`.
"""
function countyz(pstr::PauliString)
    return countyz(pstr.operator)
end

"""
    countyz(pstr::PauliStringType)

Function to count the number of Y and Z Paulis in an integer Pauli string.
"""
function countyz(pstr::PauliStringType)
    return _countbityz(pstr)
end

"""
    countyz(psum::PauliSum)

Function to count the number of Y and Z Paulis in a `PauliSum`. Returns an array of counts.
"""
function countyz(psum::PauliSum)
    return countxy(psum.op_dict)
end

"""
    countyz(psum::Dict)

Function to count the number of Y and Z Paulis in a Pauli sum in the form of a `Dict`. Returns an array of counts.
"""
function countyz(psum::Dict)
    return [countxy(pstr) for pstr in keys(psum)]
end

"""
    containsXorY(pstr::PauliString)

Check if a Pauli string contains an X or Y operator.
"""
containsXorY(pstr::PauliString) = containsXorY(pstr.operator)

"""
    containsXorY(pstr::PauliStringType)

Check if an integer Pauli string contains an X or Y operator.
"""
containsXorY(pstr::PauliStringType) = countxy(pstr) > 0

"""
    containsXorY(pstr::PauliString)

Check if a Pauli string contains a Y or Z operator.
"""
containsYorZ(pstr::PauliString) = containsYorZ(pstr.operator)

"""
    containsYorZ(pstr::PauliStringType)

Check if an integer Pauli string contains a Y or Z operator.
"""
containsYorZ(pstr::PauliStringType) = countyz(pstr) > 0


### All the commutation check functions
"""
    commutes(pstr1::PauliString, pstr2::PauliString)

Check if two Pauli strings of type `PauliString` commute.
"""
function commutes(pstr1::PauliString, pstr2::PauliString)
    return commutes(pstr1.operator, pstr2.operator)
end

"""
    commutes(pstr1::PauliStringType, pstr2::PauliStringType)

Check if two integer Pauli strings commute.
"""
function commutes(pstr1::PauliStringType, pstr2::PauliStringType)
    return _bitcommutes(pstr1, pstr2)
end

"""
    commutes(psum1::PauliSum, psum2::PauliSum)

Check if two Pauli sums of type `PauliSum` commute.
"""
function commutes(psum1::PauliSum, psum2::PauliSum)
    comm = commutator(psum1.op_dict, psum2.op_dict)
    return isempty(comm)
end

"""
    function commutes(psum1::Dict{OpType,CoeffType1}, psum2::Dict{OpType,CoeffType2}) where {OpType<:PauliStringType,CoeffType1,CoeffType2}

Check if two Pauli sums of type `PauliSum` commute.
"""
function commutes(psum1::Dict{OpType,CoeffType1}, psum2::Dict{OpType,CoeffType2}) where {OpType<:PauliStringType,CoeffType1,CoeffType2}
    comm = commutator(psum1, psum2)
    return isempty(comm)
end

"""
    commutes(pauli1::Symbol, pauli2::PauliType)

Check if two Paulis commute where one is a `Symbol` and the other is in the integer representation.
"""
function commutes(pauli1::Symbol, pauli2::PauliType)
    return commutes(pauli1, inttosymbol(pauli2))
end

"""
    commutes(pauli1::PauliType, pauli2::Symbol)

Check if two Paulis commute where one is in the integer representation and the other is a `Symbol`.
"""
function commutes(pauli1::PauliType, pauli2::Symbol)
    return commutes(pauli2, pauli1)
end

"""
    commutes(pauli1::Symbol, pauli2::Symbol)

Check if two Paulis of type `Symbol` commute.
"""
function commutes(pauli1::Symbol, pauli2::Symbol)
    if pauli1 == :I || pauli2 == :I
        return true
    else
        return pauli1 == pauli2
    end
end

## Commutator
"""
    commutator(psum1::PauliSum, psum2::PauliSum)
    
Calculate the commutator of two `PauliSum`s.
"""
function commutator(psum1::PauliSum, psum2::PauliSum)
    new_op_dict = commutator(psum1.op_dict, psum2.op_dict)
    return PauliSum(psum1.nqubits, new_op_dict)
end

"""
    commutator(pstr1::PauliString, pstr2::PauliString)

Calculate the commutator of two `PauliString`s.
"""
function commutator(pstr1::PauliString, pstr2::PauliString)
    new_coeff, new_op = commutator(pstr1.operator, pstr2.operator)
    return PauliString(pstr1.nqubits, new_op, new_coeff)
end

"""
    commutator(psum::PauliSum, pstr::PauliString)

Calculate the commutator of a `PauliSum` and a `PauliString`.
"""
commutator(psum::PauliSum, pstr::PauliString) = commutator(psum, PauliSum(pstr))

"""
    commutator(pstr::PauliString, psum::PauliSum)

Calculate the commutator of a `PauliString` and a `PauliSum`.
"""
commutator(pstr::PauliString, psum::PauliSum) = commutator(PauliSum(pstr), psum)

"""
    commutator(pstr1::PauliStringType, pstr2::PauliStringType)

Calculate the commutator of two integer Pauli strings.
Returns a tuple of the coefficient and the potentially integer Pauli string.
The coefficient is zero if the Pauli strings commute.
"""
function commutator(pstr1::PauliStringType, pstr2::PauliStringType)
    # TODO: adapt order of outputs.
    if commutes(pstr1, pstr2)
        total_sign = ComplexF64(0.0)
        new_oper = zero(typeof(pstr1))
    else
        total_sign, new_oper = pauliprod(pstr1, pstr2)
    end
    # commutator is [A, B] = AB - BA = 2AB for non-commuting (meaning anti-commuting) Paulis
    return 2 * total_sign, new_oper
end


"""
    commutator(psum1::Dict{OpType,CoeffType1}, psum2::Dict{OpType,CoeffType2}) where {OpType<:PauliStringType,CoeffType1,CoeffType2}

Calculate the commutator of two Pauli sums in the form of a `Dict`.
"""
function commutator(psum1::Dict{OpType,CoeffType1}, psum2::Dict{OpType,CoeffType2}) where {OpType<:PauliStringType,CoeffType1,CoeffType2}
    # different types of coefficients are allowed but not different types of operators

    new_pauli_dict = Dict{OpType,ComplexF64}()

    for (pauli1, coeff1) in psum1, (pauli2, coeff2) in psum2
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
"""
    pauliprod(pstr1::PauliString, pstr2::PauliString)

Calculate the product of two `PauliString`s. For example `X*Y = iZ`.
"""
function pauliprod(pstr1::PauliString, pstr2::PauliString)
    _checknumberofqubits(pstr1, pstr2)
    sign, coeff = pauliprod(pstr1.operator, pstr2.operator)
    return PauliString(pstr1.nqubits, coeff, sign * pstr1.coeff * pstr2.coeff)
end

"""
    pauliprod(pauli1::PauliStringType, pauli2::PauliStringType)

Calculate the product of two integer Pauli strings.
"""
function pauliprod(pauli1::PauliStringType, pauli2::PauliStringType)
    # This function is for when we need to globally check the sign of the product (like in general products of Paulis, not local Pauli gates)
    pauli3 = _bitpaulimultiply(pauli1, pauli2)
    sign = calculatesign(pauli1, pauli2, pauli3)
    return sign, pauli3
end

"""
    pauliprod(pauli1::Symbol, pauli2::PauliType)

Calculate the product of two Paulis where one is a `Symbol` and the other is in the integer representation.
"""
function pauliprod(pauli1::Symbol, pauli2::PauliType)
    # assume that just one qubit is involved because we check commutation with a single Symbol
    return pauliprod(symboltoint(pauli1), pauli2, 1:1)
end

"""
    pauliprod(pauli1::PauliType, pauli2::Symbol)

Calculate the product of two Paulis where one is in the integer representation and the other is a `Symbol`.
"""
function pauliprod(pauli1::PauliType, pauli2::Symbol)
    return pauliprod(pauli2, pauli1)
end

"""
    pauliprod(pauli1::Symbol, pauli2::Symbol)

Calculate the product of two Paulis of type `Symbol`.
"""
function pauliprod(pauli1::Symbol, pauli2::Symbol)
    # assume that just one qubit is involved because we check commutation with a single Symbol
    return pauliprod(symboltoint(pauli1), symboltoint(pauli2), 1:1)
end

"""
    pauliprod(pauli1::PauliType, pauli2::PauliType)

Calculate the product of two integer Paulis. 
Indicate via `changed_indices` which qubit sites to check. It can be any iterable.
"""
function pauliprod(pauli1::PauliStringType, pauli2::PauliStringType, changed_indices)
    # Calculate the Pauli product when you know on which sites the Paulis differ (changed_indices)
    pauli3 = _bitpaulimultiply(pauli1, pauli2)
    sign = calculatesign(pauli1, pauli2, pauli3, changed_indices)
    return sign, pauli3
end

"""
    calculatesign(pauli1::PauliStringType, pauli2::PauliStringType)

Calculate the sign of the product of two integer Pauli strings. Outcomes are either ±1 or ±i.
"""
function calculatesign(pauli1::PauliStringType, pauli2::PauliStringType)
    return calculatesign(pauli1, pauli2, _bitpaulimultiply(pauli1, pauli2))
end

"""
    calculatesign(pauli1::PauliStringType, pauli2::PauliStringType, pauli3::PauliStringType)

Calculate the sign of the product of two integer Pauli strings. Outcomes are either ±1 or ±i.
Takes the product of the Paulis `pauli3` as argument for efficiency. 
"""
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

"""
    calculatesign(pauli1::PauliStringType, pauli2::PauliStringType, pauli3::PauliStringType, changed_indices)

Calculate the sign of the product of two integer Pauli strings. Outcomes are either ±1 or ±i.
Takes the product of the Paulis as argument for efficiency. 
Indicate via `changed_indices` which qubit sites to check. It can be any iterable.
"""
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

"""
    generalizedlevicivita(pauli1::PauliType, pauli2::PauliType, pauli3::PauliType)

Calculate the sign of the product of two integer Paulis. Outcomes are either ±1 or ±i.
Takes the product of the Paulis as argument for efficiency. 
Indicate via `changed_indices` which qubit sites to check. It can be any iterable.

Note, this function is the foundation of `calculatesign` but assumes that the only (potentially) non-identity Pauli is on the first site.
"""
function generalizedlevicivita(pauli1::PauliType, pauli2::PauliType, pauli3::PauliType)
    # acts like levicivita but yields the correct sign for products with I or P^2, and takes care of the imaginary coefficients in Pauli products
    return generalized_levicivita_matrix[pauli1+1, pauli2+1, pauli3+1]
end

const generalized_levicivita_matrix = permutedims(cat(
        [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1], # first arg is I
        [0 1 0 0; 1 0 0 0; 0 0 0 1im; 0 0 -1im 0], # first arg is X
        [0 0 1 0; 0 0 0 -1im; 1 0 0 0; 0 1im 0 0], # first arg is Y
        [0 0 0 1; 0 0 1im 0; 0 -1im 0 0; 1 0 0 0]; # first arg is Z
        dims=3), (2, 3, 1))