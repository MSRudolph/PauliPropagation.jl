import Base: *
import Base: /
import Base: +
import Base: -
import Base: ==
import Base.length
import Base.show
import Base.sizehint!

"""
    PauliString(nqubits::Int, term::TermType, coeff::CoeffType)

`PauliString` is a `struct` that represents a Pauli string on `nqubits` qubits.
Commonly `term` is an unsigned Integer. 
See the other `PauliString` constructors for higher-level usage. 
"""
struct PauliString{TT<:PauliStringType,CT}
    nqubits::Int
    term::TT
    coeff::CT

    function PauliString(nqubits::Int, term::TT, coeff::CT) where {TT<:PauliStringType,CT}
        if nqubits <= 0
            throw(ArgumentError("Number of qubits must be a positive integer."))
        end

        return new{TT,CT}(nqubits, term, coeff)
    end
end

"""
    PauliString(nqubits::Int, pauli::Symbol, qind::Integer, coeff=1.0)
    PauliString(nqubits::Int, paulis::Vector{Symbol}, qinds::Vector{Integer}, coeff=1.0)

Constructor for a `PauliString` on `nqubits` qubits from a `Symbol` (:I, :X, :Y, :Z) or `Vector{Symbol}`.
Provide the index or indices for those symbols as `qind` or `qinds`.
The coefficient of the Pauli string in the Pauli sum defaults to 1.0.
"""
function PauliString(nqubits::Int, paulis, qinds, coeff=1.0)

    # `symboltoint()` also does checks, but we want to catch them here.
    term = symboltoint(nqubits, paulis, qinds)

    # In case an integer is passed, we probably want it to be a float.
    coeff = float(coeff)

    return PauliString(nqubits, term, coeff)
end

"""
    paulitype(pstr::PauliString)

Get the Pauli integer type of a `PauliString`.
"""
function paulitype(pstr::PauliString)
    return typeof(pstr.term)
end

"""
    coefftype(pstr::PauliString)

Get the coefficient type of a `PauliString`.
"""
function coefftype(pstr::PauliString)
    return typeof(pstr.coeff)
end

"""
    numcoefftype(pstr::PauliString)

Get the type of the numerical coefficient of a `PauliString`. 
Will return the type of the output of  `tonumber(pstr.coeff)`.
"""
function numcoefftype(pstr::PauliString)
    return typeof(tonumber(pstr.coeff))
end

"""
    *(pstr::PauliString, c::Number)

Multiply a `PauliString` by a scalar `c`. Returns a new `PauliString`.
"""
*(pstr::PauliString, c::Number) = PauliString(pstr.nqubits, pstr.term, pstr.coeff * c)
*(c::Number, pstr::PauliString) = pstr * c

"""
    /(pstr::PauliString, c::Number)

Divide a `PauliString` by a scalar `c`. Returns a new `PauliString`.
"""
/(pstr::PauliString, c::Number) = pstr * (1 / c)


# Pretty print for `PauliString`.
function show(io::IO, pstr::PauliString)
    pauli_string = inttostring(pstr.term, pstr.nqubits)
    if length(pauli_string) > 20
        pauli_string = pauli_string[1:20] * "..."
    end
    if isa(pstr.coeff, Number)
        coeff_str = round(pstr.coeff, sigdigits=5)
    elseif isa(pstr.coeff, PathProperties) && hasfield(typeof(pstr.coeff), :coeff)
        PProp = string(typeof(pstr.coeff).name.name)
        if isa(pstr.coeff.coeff, Number)
            coeff_str = "$PProp($(round(pstr.coeff.coeff, sigdigits=5)))"
        else
            coeff_str = "$PProp($(typeof(pstr.coeff.coeff)))"
        end
    else
        coeff_str = "$(typeof(pstr.coeff))"
    end
    print(io, "PauliString(nqubits: $(pstr.nqubits), $(coeff_str) * $(pauli_string))")
end


"""
    PauliSum(nqubits::Int, terms::Dict)

`PauliSum` is a `struct` that represents a sum of Pauli strings acting on `nqubits` qubits.
It is a wrapper around a dictionary Dict(Pauli string => coefficient}, where the Pauli strings are typically unsigned Integers for efficiency reasons.
"""
struct PauliSum{TT,CT}
    nqubits::Int
    terms::Dict{TT,CT}

    function PauliSum(nq::Int, terms::Dict{TT,CT}) where {TT,CT}
        if nq <= 0
            throw(ArgumentError("Number of qubits must be a positive integer."))
        end

        return new{TT,CT}(nq, terms)
    end
end

"""
    PauliSum(nqubits::Integer)

Contructor for an empty `PauliSum` on `nqubits` qubits. Element type defaults for Float64.
"""
PauliSum(nqubits::Int) = PauliSum(Float64, nqubits)

"""
    PauliSum(CoeffType, nq::Int)

Contructor for an empty `PauliSum` on `nqubits` qubits. The type of the coefficients can be provided.
"""
function PauliSum(::Type{CT}, nq::Int) where {CT}
    TT = getinttype(nq)
    return PauliSum(nq, Dict{TT,CT}())
end


"""
    PauliSum(pstr::PauliString)

Constructor for a `PauliSum` on `nqubits` qubits from a `PauliString`.
"""
PauliSum(pstr::PauliString) = PauliSum(pstr.nqubits, pstr)

"""
    PauliSum(nq::Integer, pstr::PauliString)

Constructor for a `PauliSum` on `nqubits` qubits from a `PauliString`.
"""
function PauliSum(nq::Integer, pstr::PauliString{TT,CT}) where {TT,CT}
    _checknumberofqubits(nq, pstr)
    return PauliSum(nq, Dict{TT,CT}(pstr.term => pstr.coeff))
end

"""
    PauliSum(pstrs::Vector{PauliString})

Constructor for a `PauliSum` on `nqubits` qubits from a `PauliString`.
"""
function PauliSum(pstrs::Union{AbstractArray,Tuple,Base.Generator})
    nq = _checknumberofqubits(pstrs)

    CType = promote_type(coefftype.(pstrs)...)
    psum = PauliSum(CType, nq)

    sizehint!(psum, length(pstrs))

    for pstr in pstrs
        add!(psum, pstr)
    end
    return psum
end

"""
    paulis(psum::PauliSum)

Returns an iterator over the integer pauli strings of a `PauliSum`.
Call `topaulistrings` to receive entries as `PauliString`s.
"""
function paulis(psum::PauliSum)
    return keys(psum.terms)
end

"""
    coefficients(psum::PauliSum)

Returns an iterator over the coefficients of a `PauliSum`.
Call `topaulistrings` to receive entries as `PauliString`s.
"""
function coefficients(psum::PauliSum)
    return values(psum.terms)
end

"""
    paulitype(psum::PauliSum)

Get the Pauli integer type of a `PauliSum`.
"""
function paulitype(psum::PauliSum)
    return keytype(psum.terms)
end

"""
    coefftype(psum::PauliSum)

Get the coefficient type of a `PauliSum`.
"""
function coefftype(psum::PauliSum)
    return valtype(psum.terms)
end

"""
    numcoefftype(psum::PauliSum)

Get the type of the numerical coefficient of a `PauliSum` by calling `numcoefftype()` on the coefficients.
If the `PauliSum` is empty, an error is thrown because the type cannot be inferred.
"""
function numcoefftype(psum::PauliSum)
    if length(psum) == 0
        throw(
            "Numeric coefficient type cannot be inferred from an empty PauliSum." *
            "Consider defining a `numcoefftype(psum::$(typeof(psum)))` method.")
    end
    return typeof(tonumber(first(coefficients(psum))))
end


"""
    getcoeff(psum::PauliSum, pstr::Integer)

Get the coefficient of an integer Pauli string in a `PauliSum`. Defaults to 0 if the Pauli string is not in the `PauliSum`.
Requires that the integer Pauli string `pstr` is the same type as the integer Pauli strings in `psum`.
"""
function getcoeff(psum::PauliSum{TT1,CT}, pstr::TT2) where {TT1,TT2,CT}
    return get(psum.terms, pstr, zero(numcoefftype(psum)))
end

"""
    getcoeff(psum::PauliSum, pstr::PauliString)

Get the coefficient of a `PauliString` in a `PauliSum`. Defaults to 0 if the Pauli string is not in the `PauliSum`.
Requires that the integer Pauli string in `pstr` is the same type as the integer Pauli strings in `psum`.
"""
function getcoeff(psum::PauliSum{TT1,CT1}, pstr::PauliString{TT2,CT2}) where {TT1,TT2,CT1,CT2}
    return getcoeff(psum, pstr.term)
end


"""
    getcoeff(psum::PauliSum, pauli::Symbol, qind::Integer)

Get the coefficient of a Pauli string in a `PauliSum` by providing the Pauli string as a Symbol acting on qubit `qind`. 
This is consistent with how Pauli strings can be added to a `PauliSum` via `add!()`. 
Defaults to 0 if the Pauli string is not in the `PauliSum`.
"""
function getcoeff(psum::PauliSum, pauli::Symbol, qind::Integer)
    return getcoeff(psum, symboltoint(psum.nqubits, pauli, qind))
end

"""
    getcoeff(psum::PauliSum, pstr::Vector{Symbol}, qinds::Vector{Int})

Get the coefficient of a Pauli string in a `PauliSum` by providing the Pauli string `pstr` as a vector of Symbols acting on qubits `qinds`. 
This is consistent with how Pauli strings can be added to a `PauliSum` via `add!()`. 
Defaults to 0 if the Pauli string is not in the `PauliSum`.
"""
function getcoeff(psum::PauliSum, pstr, qinds)
    return getcoeff(psum, symboltoint(psum.nqubits, pstr, qinds))
end

"""
    getcoeff(psum::PauliSum, pstr::Vector{Symbol})

Get the coefficient of a Pauli string in a `PauliSum` by providing the Pauli string `pstr` as a vector of Symbols acting on all qubits. 
This is consistent with how Pauli strings can be added to a `PauliSum` via `add!()`. 
Defaults to 0 if the Pauli string is not in the `PauliSum`.
"""
function getcoeff(psum::PauliSum, pstr::Vector{Symbol})
    return getcoeff(psum, symboltoint(pstr))
end

"""
    tonumber(val::Number)

Trivial function returning a numerical value of a number.
Will be overloaded for custom wrapper types like `PathProperties`.
"""
function tonumber(val::Number)
    return val
end

"""
    norm(psum::PauliSum, L=2)

Calculate the norm of a `PauliSum` with respect to the `L`-norm. 
Calls LinearAlgebra.norm on the coefficients of the `PauliSum`.
"""
function norm(psum::PauliSum, L=2)
    if length(psum) == 0
        return 0.0
    end
    return LinearAlgebra.norm((tonumber(coeff) for coeff in coefficients(psum)), L)
end

"""
    topaulistrings(psum::PauliSum)

Returns the Pauli strings in a `PauliSum` and their coefficients as a list of `PauliString`.
"""
topaulistrings(psum::PauliSum) = [PauliString(psum.nqubits, pauli, coeff) for (pauli, coeff) in psum.terms]

"""
Copy a `PauliSum` by copying its `terms` field.
"""
Base.copy(psum::PauliSum) = PauliSum(psum.nqubits, copy(psum.terms))

"""
Iterator for `PauliSum` returns an iterator over (pstr, coeff) pairs`.
"""
Base.iterate(psum::PauliSum, state=1) = iterate(psum.terms, state)


# Pretty print for `PauliSum`.
function show(io::IO, psum::PauliSum)
    if length(psum.terms) == 0
        dict_string = "(no Pauli strings)"
    else
        dict_string = _getprettystr(psum.terms, psum.nqubits)
    end
    print(io, "PauliSum(nqubits: $(psum.nqubits), $dict_string)")
end


"""
    length(psum::PauliSum)

Number of terms in the `PauliSum`.
"""
length(psum::PauliSum) = length(psum.terms)

"""
    sizehint!(psum::PauliSum, n)

Hint to the `PauliSum` to reserve space for `n` terms.
"""
sizehint!(psum::PauliSum, n) = sizehint!(psum.terms, n)


"""
    ==(psum1::PauliSum, psum2::PauliSum)

Equality check for `PauliSum`.
"""
function ==(psum1::PauliSum, psum2::PauliSum)
    if psum1.nqubits != psum2.nqubits
        return false
    end

    return psum1.terms == psum2.terms
end

"""
    ≈(psum1::PauliSum, psum2::PauliSum)

Approximate equality check for `PauliSum`.
Simply calls `isapprox()` on the coefficients of the contained Pauli strings.
"""
function Base.:≈(psum1::PauliSum{TT1,CT1}, psum2::PauliSum{TT2,CT2}) where {TT1,CT1,TT2,CT2}
    if TT1 != TT2
        return false
    end

    if psum1.nqubits != psum2.nqubits
        return false
    end

    # we don't strictly need to check the length of the dictionaries
    # small values are allowed for Pauli strings that don't exist in both
    # our solution is to check for approximate equality both ways
    for (pstr, coeff) in psum1
        if !isapprox(get(psum2.terms, pstr, CT2(0.0)), coeff)
            return false
        end
    end
    for (pstr, coeff) in psum2
        if !isapprox(get(psum1.terms, pstr, CT1(0.0)), coeff)
            return false
        end
    end

    return true
end


## Arithmetic operations
# They all deepcopy

"""
    *(psum::PauliSum, c::Number)

Multiply a `PauliSum` by a scalar `c`. This copies the PauliSum.
"""
function *(psum::PauliSum, c::Number)
    ps_copy = deepcopy(psum)
    mult!(ps_copy, c)
    return ps_copy
end

*(c::Number, psum::PauliSum) = psum * c


"""
    /(psum::PauliSum, c::Number)

Divide a `PauliSum` by a scalar `c`. This copies the PauliSum.
"""
function /(psum::PauliSum, c::Number)
    ps_copy = deepcopy(psum)
    mult!(ps_copy, 1 / c)
    return ps_copy
end

"""
    +(psum::PauliSum, c::Number)
    +(c::Number, psum::PauliSum)

Addition of c * Identity to a `PauliSum`. This copies the PauliSum.
"""
function +(psum::PauliSum{TT,CT}, c::Number) where {TT,CT}
    psum = deepcopy(psum)
    add!(psum, identitypauli(TT), c)
    return psum
end

+(c::Number, psum::PauliSum) = psum + c


"""
    +(pstr1::PauliString, pstr2::PauliString)

Addition of two `PauliString`s. Returns a PauliSum.
"""
function +(pstr1::PauliString{TT,CT1}, pstr2::PauliString{TT,CT2}) where {TT,CT1,CT2}
    nq = _checknumberofqubits(pstr1, pstr2)

    # get a compatibel coefficient type
    CType = promote_type(coefftype(pstr1), coefftype(pstr2))
    psum = PauliSum(CType, nq)
    add!(psum, pstr1)
    add!(psum, pstr2)
    return psum
end

"""
    +(pstr::PauliString, psum::PauliSum)
    +(psum::PauliSum, pstr::PauliString)

Addition of a `PauliString` to a `PauliSum`. Returns a `PauliSum`.
"""
function +(psum::PauliSum{TT,CT1}, pstr::PauliString{TT,CT2}) where {TT,CT1,CT2}
    nq = _checknumberofqubits(psum, pstr)

    # get a compatible coefficient type
    CType = promote_type(CT1, CT2)
    new_psum = PauliSum(CType, nq)

    add!(new_psum, psum)
    add!(new_psum, pstr)
    return new_psum
end

+(pstr::PauliString, psum::PauliSum) = psum + pstr


"""
    +(psum1::PauliSum, psum2::PauliSum)
Addition of two `PauliSum`s. Returns a `PauliSum`.
"""
function +(psum1::PauliSum{TT1,CT1}, psum2::PauliSum{TT2,CT2}) where {TT1,TT2,CT1,CT2}

    # throw custom error if paulitypes are not the same
    _checktermtype(psum1, psum2)
    nq = _checknumberofqubits(psum1, psum2)

    # get a compatible coefficient type
    CType = promote_type(CT1, CT2)
    psum = PauliSum(CType, nq)

    add!(psum, psum1)
    add!(psum, psum2)
    return psum
end



"""
    -(pstr1::PauliString, pstr2::PauliString)

Subtract two `PauliString`s. Returns a PauliSum.
"""
function -(pstr1::PauliString, pstr2::PauliString)
    return pstr1 + (-1 * pstr2)
end

"""
    -(pstr::PauliString, psum::PauliSum)
    -(psum::PauliSum, pstr::PauliString)

Subtract a `PauliString` from a `PauliSum` or vice versa.
Returns a `PauliSum`.
"""
function -(psum::PauliSum, pstr::PauliString)
    return psum + (-1 * pstr)
end

-(pstr::PauliString, psum::PauliSum) = psum - pstr


"""
    -(psum1::PauliSum, psum2::PauliSum)

Subtract two `PauliSum`s. Returns a `PauliSum`.
"""
function -(psum1::PauliSum, psum2::PauliSum)
    deepcopy(psum2)
    mult!(psum2, -1)
    return psum1 + psum2
end


# Pauli products
"""
    *(pstr1::PauliString, pstr2::PauliString)

Perform a Pauli product of two `PauliString`s. 
"""
function *(pstr1::PauliString, pstr2::PauliString)
    _checktermtype(pstr1, pstr2)
    _checknumberofqubits(pstr1, pstr2)

    return pauliprod(pstr1, pstr2)
end

"""
    *(pstr::PauliString, psum::PauliSum)
    *(psum::PauliSum, pstr::PauliString)

Perform a Pauli product of a `PauliString` with a `PauliSum`.
Returns a `PauliSum` with complex coefficients.
"""
function *(psum::PauliSum, pstr::PauliString)
    _checktermtype(psum, pstr)
    nq = _checknumberofqubits(psum, pstr)

    new_psum = PauliSum(ComplexF64, nq)
    sizehint!(new_psum, length(psum))

    for (term, coeff) in psum
        new_pstr = pauliprod(PauliString(nq, term, coeff), pstr)
        add!(new_psum, new_pstr)
    end
    return new_psum
end

function *(pstr::PauliString, psum::PauliSum)
    _checktermtype(psum, pstr)
    nq = _checknumberofqubits(psum, pstr)

    # TODO: this is literally the same, just argument oder reversed in pauliprod()
    new_psum = PauliSum(ComplexF64, nq)
    sizehint!(new_psum, length(psum))

    for (term, coeff) in psum
        new_pstr = pauliprod(pstr, PauliString(nq, term, coeff))
        add!(new_psum, new_pstr)
    end
    return new_psum
end


"""
    *(psum1::PauliSum, psum2::PauliSum)

Perform a Pauli product of two `PauliSum`s.
Returns a `PauliSum` with complex coefficients.
"""
function *(psum1::PauliSum, psum2::PauliSum)
    _checktermtype(psum1, psum2)
    nq = _checknumberofqubits(psum1, psum2)

    psum = PauliSum(ComplexF64, nq)
    sizehint!(psum, length(psum1))

    for (pstr1, coeff1) in psum1
        for (pstr2, coeff2) in psum2
            pstr, sign = pauliprod(pstr1, pstr2)
            coeff = coeff1 * coeff2 * sign
            add!(psum, pstr, coeff)
        end
    end
    return psum

end


## In-place Multiplication

"""
    mult!(psum::PauliSum, c::Number)

Multiply a `PauliSum` by a scalar `c` in-place.
"""
function mult!(psum::PauliSum, c::Number)
    # multiply in-place
    for (k, v) in psum.terms
        psum.terms[k] *= c
    end
    return psum
end

## In-place Addition

"""
    add!(psum::PauliSum, pauli::Symbol, qind::Integer, coeff=1.0)
    add!(psum::PauliSum, paulis::Vector{Symbol}, qinds::Vector{Integer}, coeff=1.0)

Add a Pauli string to a `PauliSum` `psum`. Changes `psum` in-place.
Provide the Pauli string as a `Symbol` (:I, :X, :Y, :Z) or `Vector{Symbol}`.
Provide the index or indices for those symbols as `qind` or `qinds`.
The coefficient of the Pauli string in the Pauli sum defaults to 1.0.
"""
function add!(psum::PauliSum, paulis::Union{Symbol,Vector{Symbol}}, qinds, coeff=coefftype(psum)(1.0))
    return add!(psum, PauliString(psum.nqubits, paulis, qinds, coeff))
end

"""
    add!(psum::PauliSum, pstr::PauliString)

Add a `PauliString` `pstr` to a `PauliSum` `psum`. Changes `psum` in-place.
`psum` and `pstr` need to be defined on the same number of qubits and have the same coefficient type.
"""
function add!(psum::PauliSum{TT1,CT1}, pstr::PauliString{TT2,CT2}) where {TT1,TT2,CT1,CT2}
    _checktermtype(psum, pstr)
    _checknumberofqubits(psum, pstr)

    # this is supposed to error if pstr.coeff cannot be converted to CT1
    # because this is an in-place operation
    pstr_coeff = convert(CT1, pstr.coeff)
    add!(psum, pstr.term, pstr_coeff)
    return psum
end

"""
    add!(psum1::PauliSum, psum2::PauliSum)

Add two `PauliSum`s `psum1` and `psum2`. Changes `psum1` in-place.
`psum1` and `psum2` need to be defined on the same number of qubits and have the same coefficient type.
"""
function add!(psum1::PauliSum{TT1,CT1}, psum2::PauliSum{TT2,CT2}) where {TT1,TT2,CT1,CT2}
    _checktermtype(psum1, psum2)
    _checknumberofqubits(psum1, psum2)

    add!(psum1.terms, psum2.terms)
    return psum1
end


"""
    add!(psum::PauliSum{Integer, CoeffType}, pstr::Integer, coeff::CoeffType)

Add a Pauli string `pstr` with coefficient `coeff` to a `PauliSum` `psum`. This changes `psum` in-place.
`pstr` needs to have the same type as `paulitype(psum)`, and `coeff` needs to have the same type as `coefftype(psum)`.
"""
function add!(psum::PauliSum{TT,CT1}, pstr::TT, coeff::CT2) where {TT,CT1,CT2}
    add!(psum.terms, pstr, coeff)
    return psum
end

function add!(psum1::Dict{TT,CT1}, psum2::Dict{TT,CT2}) where {TT,CT1,CT2}
    ## Lower level addition of two Pauli sum dictionaries
    for (pstr, coeff) in psum2
        add!(psum1, pstr, coeff)
    end
    return psum1
end


function add!(psum::Dict{TT,CT1}, pstr::TT, coeff::CT2) where {TT,CT1,CT2}
    ## Lower level addition of a Pauli string to a Pauli sum dictionary

    # don't add if the coefficient is 0
    if tonumber(coeff) == 0
        return psum
    end

    if haskey(psum, pstr)
        new_coeff = psum[pstr] + coeff
        if tonumber(new_coeff) == 0
            delete!(psum, pstr)
        else
            psum[pstr] = new_coeff
        end

    else
        psum[pstr] = coeff
    end
    return psum
end


## Set in Pauli sum
"""
    set!(psum::PauliSum{TermType, CoeffType}, pstr::TermType, coeff::CoeffType)

In-place setting the coefficient of a Pauli string in a `PauliSum` dictionary.
The type of the Pauli string needs to be the keytype=`TermType` of the dictionary, and the coefficient `coeff` needs to be the valuetype=`CoeffType`.
"""
function set!(psum::PauliSum{TT,CT1}, pstr::TT, coeff::CT2) where {TT,CT1,CT2}
    set!(psum.terms, pstr, coeff)
    return psum
end


function set!(psum::Dict{TT,CT1}, pstr::TT, coeff::CT2) where {TT,CT1,CT2}
    # lower-level set!() for Pauli sum dict

    # delete if the coefficient would be set to 0
    if tonumber(coeff) == 0
        delete!(psum, pstr)

    else
        psum[pstr] = coeff
    end
    return psum
end

## Helper functions
function Base.delete!(psum::PauliSum{TT,CT}, pstr::TT) where {TT,CT}
    delete!(psum.terms, pstr)
    return psum
end

"""
    empty!(psum::PauliSum)

Empty the `PauliSum` by emptying the dictionary on the `terms` fields. 
"""
Base.empty!(psum::PauliSum) = empty!(psum.terms)


"""
    similar(psum::PauliSum)

Create a new `PauliSum` with the same number of qubits and coefficient type as `psum`.
Calls `sizehint!()` with `length(psum)` on the dictionary of the new `PauliSum`. 
"""
function Base.similar(psum::PauliSum)
    return PauliSum(psum.nqubits, similar(psum.terms))
end

function Base.similar(psum::Dict{TT,CT}) where {TT,CT}
    new_psum = Dict{TT,CT}()
    sizehint!(new_psum, length(psum))
    return new_psum
end



# Checks whether the number of qubits `nqubits` is the same between our datatypes.
function _checknumberofqubits(nqubits::Int, pobj::Union{PauliString,PauliSum})
    if nqubits != pobj.nqubits
        throw(
            ArgumentError(
                "Number of qubits ($(nqubits)) must equal number of qubits ($(pobj.nqubits)) in $(typeof(pobj))"
            )
        )
    end
    return nqubits
end


# Checks whether the number of qubits `nqubits` is the same between as the length of the vector `pstr`.
function _checknumberofqubits(nqubits::Int, pstr)
    if nqubits != length(pstr)
        throw(
            ArgumentError(
                "Number of qubits ($(op1.nqubits)) must equal number of qubits ($(length(pstr))) in $(typeof(pstr))"
            )
        )
    end
    return nqubits
end


# Checks whether the number of qubits `nqubits` is the same between our datatypes.
function _checknumberofqubits(pobj1::Union{PauliString,PauliSum}, pobj2::Union{PauliString,PauliSum})
    if pobj1.nqubits != pobj2.nqubits
        throw(
            ArgumentError(
                "Number of qubits ($(pobj1.nqubits)) in $(typeof(pobj1)) must equal number of qubits ($(pobj2.nqubits)) in $(typeof(pobj2))"
            )
        )
    end
    return pobj1.nqubits
end

"""
Checks whether the number of qubits `nqubits` is the same between in some collection.
"""
function _checknumberofqubits(pobjects::Union{AbstractArray,Tuple,Base.Generator})

    if !allequal(pobj.nqubits for pobj in pobjects)
        throw(
            ArgumentError(
                "Number of qubits in passed collection of type $(typeof(pobjects)) is not consistent."
            )
        )
    end
    return first(pobjects).nqubits
end


# throw error for miss-matched term/Pauli types

function _checktermtype(pobj1, pobj2)
    if paulitype(pobj1) != paulitype(pobj2)
        throw(ArgumentError("Pauli types do not match. Got $(TT1) and $(TT2)."))
    end
end