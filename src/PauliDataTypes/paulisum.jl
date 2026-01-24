

"""
    PauliSum(nqubits::Int, terms::Dict) <: AbstractPauliSum

`PauliSum` is a `struct` that represents a sum of Pauli strings acting on `nqubits` qubits.
It is a wrapper around a dictionary Dict(Pauli string => coefficient}, where the Pauli strings are typically unsigned Integers for efficiency reasons.
"""
struct PauliSum{TT,CT} <: AbstractPauliSum
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

PropagationBase.storage(psum::PauliSum) = psum.terms

"""
    nqubits(psum::PauliSum)

Get the number of qubits that the `PauliSum` is defined on.
"""
nqubits(psum::PauliSum) = psum.nqubits


"""
Copy a `PauliSum` by copying its `terms` field.
"""
Base.copy(psum::PauliSum) = PauliSum(psum.nqubits, copy(psum.terms))


# Pretty print for `PauliSum`.
function Base.show(io::IO, psum::PauliSum)
    if length(psum.terms) == 0
        dict_string = "(no Pauli strings)"
    else
        dict_string = _getprettystr(psum.terms, psum.nqubits)
    end
    print(io, "PauliSum(nqubits: $(psum.nqubits), $dict_string)")
end


"""
    sizehint!(psum::PauliSum, n)

Hint to the `PauliSum` to reserve space for `n` terms.
"""
Base.sizehint!(psum::PauliSum, n) = sizehint!(psum.terms, n)


# TODO: Move these to PropagationBase
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

-(pstr::PauliString, psum::PauliSum) = mult!(psum - pstr, -1)


"""
    -(psum1::PauliSum, psum2::PauliSum)

Subtract two `PauliSum`s. Returns a `PauliSum`.
"""
function -(psum1::PauliSum, psum2::PauliSum)
    psum2 = deepcopy(psum2)
    mult!(psum2, -1)
    return psum1 + psum2
end


# Pauli products

"""
    *(pstr::PauliString, psum::PauliSum)
    *(psum::PauliSum, pstr::PauliString)

Perform a Pauli product of a `PauliString` with a `PauliSum`.
Returns a `PauliSum` with complex coefficients.
"""
function *(psum::PauliSum, pstr::PauliString)

    psum2 = PauliSum(pstr)
    return pauliprod(psum, psum2)
end

function *(pstr::PauliString, psum::PauliSum)

    psum2 = PauliSum(pstr)
    return pauliprod(psum2, psum)
end


"""
    *(psum1::PauliSum, psum2::PauliSum)

Perform a Pauli product of two `PauliSum`s.
Returns a `PauliSum` with complex coefficients.
"""
function *(psum1::PauliSum, psum2::PauliSum)

    return pauliprod(psum1, psum2)
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

