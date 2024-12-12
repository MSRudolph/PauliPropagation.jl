import Base: *
import Base: /
import Base: +
import Base: -

"""
    PauliString(nqubits::Int, term::TermType, coeff::CoeffType)

`PauliString` is a struct that represents a Pauli string on `nqubits` qubits.
Commonly `term` is an unsigned Integer. See the other `PauliString` constructors for details. 
"""
struct PauliString{TermType<:PauliStringType,CoeffType}
    nqubits::Int
    term::TermType
    coeff::CoeffType
end

"""
    PauliString(nqubits::Integer, pauli::Symbol, qind::Integer, coeff=1.0)

Constructor for a `PauliString` on `nqubits` qubits from a Symbol (:X, :Y, :Z) representing a single non-identity Pauli on qubit `qind` with coefficient `coeff`.
"""
function PauliString(nqubits::Integer, pauli::Symbol, qind::Integer, coeff=1.0)
    pauli = symboltoint(nqubits, pauli, qind)
    coeff = _convertcoefficients(coeff)
    return PauliString(nqubits, pauli, coeff)
end

"""
    PauliString(nqubits, pstr, qinds, coeff=1.0)

Constructor for a `PauliString` on `nqubits` qubits from a list of Symbols representing non-identity Paulis on qubits `qinds` with coefficient `coeff`.
"""
function PauliString(nqubits, pstr, qinds, coeff=1.0)
    pauli = symboltoint(nqubits, pstr, qinds)
    coeff = _convertcoefficients(coeff)
    return PauliString(nqubits, pauli, coeff)
end

"""
    term(pstr::PauliString)

Get the lower-level representation of a `PauliString`.
This returns the `term` field of the `PauliString`. 
"""
term(pstr::PauliString) = pstr.term

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
Will get the type of the `coeff` field of a potential PathProperties type.
"""
function numcoefftype(pstr::PauliString)
    return numcoefftype(pstr.coeff)
end


import Base.show
"""
Pretty print for `PauliString`.
"""
function show(io::IO, pstr::PauliString)
    pauli_string = inttostring(pstr.term, pstr.nqubits)
    if length(pauli_string) > 20
        pauli_string = pauli_string[1:20] * "..."
    end
    if isa(pstr.coeff, Number)
        coeff_str = round(pstr.coeff, sigdigits=5)
    elseif isa(pstr.coeff, PathProperties)
        if isa(pstr.coeff, Number)
            coeff_str = "PathProperty($(round(pstr.coeff, sigdigits=5)))"
        else
            coeff_str = "PathProperty($(typeof(pstr.coeff)))"
        end
    else
        coeff_str = "($(typeof(pstr.coeff)))"
    end
    print(io, "PauliString(nqubits: $(pstr.nqubits), $(coeff_str) * $(pauli_string))")
end


"""
    PauliSum(nqubits::Int, terms::Dict)

`PauliSum`` is a struct that represents a sum of Pauli strings acting on `nqubits` qubits.
It is a wrapper around a dictionary Dict(Pauli string => coefficient}, where the Pauli strings are typically unsigned Integers for efficiency reasons.
"""
struct PauliSum{TermType<:Integer,CoeffType}
    nqubits::Int
    terms::Dict{TermType,CoeffType}
end

"""
    PauliSum(nqubits::Integer)

Contructor for an empty `PauliSum` on `nqubits` qubits. Element type defaults for Float64.
"""
PauliSum(nqubits::Integer) = PauliSum(nqubits, Float64)

"""
    PauliSum(nq::Int, ELTYPE::T) where {T<:DataType}

Contructor for an empty `PauliSum` on `nqubits` qubits. Element type can be provided.
"""
function PauliSum(nq::Int, ELTYPE::T) where {T<:DataType}
    TermType = getinttype(nq)
    return PauliSum(nq, Dict{TermType,ELTYPE}())
end

"""
    PauliSum(nqubits::Integer, psum::Dict{Vector{Symbol},CoeffType}) where {CoeffType}

Constructor for a `PauliSum` on `nqubits` qubits from a dictionary of {Vector{Symbols} => coefficients}.
"""
function PauliSum(nqubits::Integer, psum::Dict{Vector{Symbol},CoeffType}) where {CoeffType}

    _checknumberofqubits.(nqubits, keys(psum))

    int_dict = Dict(symboltoint(k) => _convertcoefficients(v) for (k, v) in psum)

    return PauliSum(nqubits, int_dict)
end

"""
    PauliSum(psum::PauliSum)

Trivial constructor for a `PauliSum` on `nqubits` qubits from a `PauliSum`. Returns the same `PauliSum` and does not copy.
"""
PauliSum(psum::PauliSum) = psum

"""
    PauliSum(pstr::PauliString)

Constructor for a `PauliSum` on `nqubits` qubits from a `PauliString`.
"""
PauliSum(pstr::PauliString) = PauliSum(pstr.nqubits, pstr)

"""
    PauliSum(nq::Integer, pstr::PauliString{TermType,CoeffType}) where {TermType,CoeffType}

Constructor for a `PauliSum` on `nqubits` qubits from a `PauliString`.
"""
function PauliSum(nq::Integer, pstr::PauliString{TermType,CoeffType}) where {TermType,CoeffType}
    _checknumberofqubits(nq, pstr)
    return PauliSum(nq, Dict{TermType,CoeffType}(pstr.term => pstr.coeff))
end

"""
    terms(psum::PauliSum)

Returns the data structure of the `PauliSum` containing the Pauli strings and their coefficients.
Currently a dictionary.
"""
terms(psum::PauliSum) = psum.terms

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
    paulitype(psum::Dict)

Get the Pauli integer type of a Pauli sum dict.
"""
function paulitype(psum::Dict)
    return keytype(psum)
end

"""
    coefftype(psum::PauliSum)

Get the coefficient type of a `PauliSum`.
"""
function coefftype(psum::PauliSum)
    return valtype(psum.terms)
end

"""
    coefftype(psum::Dict)

Get the coefficient type of a `PauliSum`.
"""
function coefftype(psum::Dict)
    return valtype(psum)
end

"""
    numcoefftype(psum::PauliSum)

Get the type of the numerical coefficient of a `PauliSum`. 
Will get the type of the `coeff` field of a potential PathProperties type.
"""
function numcoefftype(psum::PauliSum)
    return numcoefftype(psum.terms)
end

"""
    numcoefftype(psum::Dict)

Get the type of the numerical coefficient of a pauli sum dict. 
Will get the type of the `coeff` field of a potential PathProperties type.
"""
function numcoefftype(psum::Dict)
    return numcoefftype(valtype(psum))
end

"""
    numcoefftype(::Number)

Get the type of the number.
"""
function numcoefftype(::T) where {T<:Number}
    return T
end

"""
    numcoefftype(::Type{Number})

Return the input type if it is a Number type.
"""
function numcoefftype(::Type{T}) where {T<:Number}
    return T
end

"""
    getcoeff(psum::PauliSum{PauliStringType,CoeffType}, pstr::PauliStringType)

Get the coefficient of an integer Pauli string in a `PauliSum`. Defaults to 0 if the Pauli string is not in the `PauliSum`.
Requires that the integer Pauli string `pstr` is the same type as the integer Pauli strings in `psum`.
"""
function getcoeff(psum::PauliSum{TermType,CoeffType}, pstr::TermType) where {TermType,CoeffType}
    # TODO: This is not yet compatible with `PathProperties`
    if CoeffType <: PathProperties
        throw(ArgumentError("This function is not yet compatible with PathProperties."))
    end
    return get(psum.terms, pstr, CoeffType(0))
end

"""
    getcoeff(psum::PauliSum, pstr::PauliString)

Get the coefficient of a `PauliString` in a `PauliSum`. Defaults to 0 if the Pauli string is not in the `PauliSum`.
Requires that the integer Pauli string in `pstr` is the same type as the integer Pauli strings in `psum`.
"""
function getcoeff(psum::PauliSum{TermType,CoeffType1}, pstr::PauliString{TermType,CoeffType2}) where {TermType,CoeffType1,CoeffType2}
    # TODO: This is not yet compatible with `PathProperties`
    if CoeffType <: PathProperties
        throw(ArgumentError("This function is not yet compatible with PathProperties."))
    end
    return get(psum.terms, pstr.term, CoeffType(0))
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


# TODO getnumcoeff() for PauliSum and PauliString


"""
    getnumcoeff(val::Number)

Trivial function returning a numerical value of a number.
"""
function getnumcoeff(val::Number)
    return val
end

# TODO: implement norm for PauliSum

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
Iterator for `PauliSum` returns the iterator over its `terms`.
"""
Base.iterate(psum::PauliSum, state=1) = iterate(psum.terms, state)


import Base.show
"""
Pretty print for `PauliSum`.
"""
function show(io::IO, psum::PauliSum)
    if length(psum.terms) == 0
        dict_string = "(no Pauli strings)"
    else
        dict_string = _getprettystr(psum.terms, psum.nqubits)
    end
    print(io, "PauliSum(nqubits: $(psum.nqubits), $dict_string)")
end

import Base.length
"""
    length(psum::PauliSum)

Number of terms in the `PauliSum`.
"""
length(psum::PauliSum) = length(psum.terms)

import Base: ==
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

"""
    *(psum::PauliSum, c::Number)

Multiply a `PauliSum` by a scalar `c`. This copies the PauliSum.
"""
function *(psum::PauliSum, c::Number)
    ps_copy = deepcopy(psum)
    mult!(ps_copy, c)
    return ps_copy
end

"""
    *(c::Number, psum::PauliSum)

Multiply a `PauliSum` by a scalar `c`. This copies the PauliSum.
"""
function *(c::Number, psum::PauliSum)
    return psum * c
end

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
    +(pstr1::PauliString{TermType,CoeffType}, pstr2::PauliString{TermType,CoeffType})

Addition of two `PauliString`s. Returns a PauliSum.
"""
function +(pstr1::PauliString{TermType,CoeffType}, pstr2::PauliString{TermType,CoeffType}) where {TermType,CoeffType}
    _checknumberofqubits(pstr1, pstr2)
    psum = PauliSum(pstr1)
    add!(psum, pstr2)
    return psum
end

"""
    +(psum::PauliSum{TermType,CoeffType}, pstr::PauliString{TermType,CoeffType})

Addition of a `PauliString` to a `PauliSum`. Returns a `PauliSum`.
"""
function +(psum::PauliSum{TermType,CoeffType}, pstr::PauliString{TermType,CoeffType}) where {TermType,CoeffType}
    _checknumberofqubits(psum, pstr)
    psum = deepcopy(psum)
    add!(psum, pstr)
    return psum
end

"""
    +(psum1::PauliSum{TermType,CoeffType}, psum2::PauliSum{TermType,CoeffType})
Addition of two `PauliSum`s. Returns a `PauliSum`.
"""
function +(psum1::PauliSum{TermType,CoeffType}, psum2::PauliSum{TermType,CoeffType}) where {TermType,CoeffType}
    _checknumberofqubits(psum1, psum2)
    psum1 = deepcopy(psum1)
    add!(psum1, psum2)
    return psum1
end

"""
    add!(psum::PauliSum{TermType,CoeffType}, pstr::PauliString{TermType,CoeffType})

Addition of a `PauliString` to a `PauliSum`. Changes the `PauliSum` in-place.
"""
function add!(psum::PauliSum{TermType,CoeffType}, pstr::PauliString{TermType,CoeffType}) where {TermType,CoeffType}
    _checknumberofqubits(psum, pstr)
    add!(psum.terms, pstr.term, pstr.coeff)
    return psum
end

"""
    add!(psum1::PauliSum{TermType,CoeffType}, psum2::PauliSum{TermType,CoeffType})

Addition of two `PauliSum`s. Changes the first `PauliSum` in-place.
"""
function add!(psum1::PauliSum{TermType,CoeffType}, psum2::PauliSum{TermType,CoeffType}) where {TermType,CoeffType}
    _checknumberofqubits(psum1, psum2)
    add!(psum1.terms, psum2.terms)
    return psum1
end

"""
    add!(psum1::Dict{TermType,CoeffType}, psum2::Dict{TermType,CoeffType})

Addition of two Pauli sum dicts. Changes the first pauli sum in-place.
"""
function add!(psum1::Dict{TermType,CoeffType}, psum2::Dict{TermType,CoeffType}) where {TermType,CoeffType}
    for (pstr, coeff) in psum2
        add!(psum1, pstr, coeff)
    end
    return psum1
end

"""
    add!(psum::PauliSum{TermType, CoeffType}, pstr::TermType, coeff::CoeffType)

In-place addition of a Pauli string `pstr` with coefficient `coeff` to a `PauliSum`.
"""
function add!(psum::PauliSum{TermType,CoeffType}, pstr::TermType, coeff::CoeffType) where {TermType,CoeffType}
    add!(psum.terms, pstr, coeff)
    return psum
end

"""
    add!(psum::Dict{TermType, CoeffType}, pstr::TermType, coeff::CoeffType)

In-place addition of a Pauli string `pstr` with coefficient `coeff` to a Pauli sum dictionary.
"""
function add!(psum::Dict{TermType,CoeffType}, pstr::TermType, coeff::CoeffType) where {TermType,CoeffType}
    # don't add if the coefficient is 0
    if getnumcoeff(coeff) == 0
        return psum
    end

    if haskey(psum, pstr)
        new_coeff = psum[pstr] + coeff
        if getnumcoeff(new_coeff) == 0
            delete!(psum, pstr)
        else
            psum[pstr] = new_coeff
        end

    else
        psum[pstr] = coeff
    end
    return psum
end

"""
    add!(psum::PauliSum, pauli::Union{Symbol, Vector{Symbol}}, qind::Union{Integer, Vector{Integer}}, coeff=1.0)

In-place addition a Pauli string to `PauliSum` by providing the Pauli string as a Symbol or a vector of symbols acting on qubits `qinds`.
Coefficient defaults to 1.0.
"""
function add!(psum::PauliSum, paulis::Union{Symbol,Vector{Symbol}}, qinds::Union{T,Vector{T}}, coeff=1.0) where {T<:Integer}
    return add!(psum, PauliString(psum.nqubits, paulis, qinds, coeff))
end

## Substraction
"""
    -(pstr1::PauliString{TermType,CoeffType}, pstr2::PauliString{TermType,CoeffType})

Subtraction of two `PauliString`s. Returns a PauliSum.
"""
function -(pstr1::PauliString{TermType,CoeffType}, pstr2::PauliString{TermType,CoeffType}) where {TermType,CoeffType}
    _checknumberofqubits(pstr1, pstr2)
    psum = PauliSum(pstr1)
    subtract!(psum, pstr2)
    return psum
end

"""
    -(psum::PauliSum{TermType,CoeffType}, pstr::PauliString{TermType,CoeffType})

Subtraction of a `PauliString` to a `PauliSum`. Returns a `PauliSum`.
"""
function -(psum::PauliSum{TermType,CoeffType}, pstr::PauliString{TermType,CoeffType}) where {TermType,CoeffType}
    _checknumberofqubits(psum, pstr)
    psum = deepcopy(psum)
    subtract!(psum, pstr)
    return psum
end

"""
    -(psum1::PauliSum{TermType,CoeffType}, psum2::PauliSum{TermType,CoeffType})

Subtract of two `PauliSum`s. Returns a `PauliSum`.
"""
function -(psum1::PauliSum{TermType,CoeffType}, psum2::PauliSum{TermType,CoeffType}) where {TermType,CoeffType}
    _checknumberofqubits(psum1, psum2)
    psum1 = deepcopy(psum1)
    subtract!(psum1, psum2)
    return psum1
end

"""
    subtract!(psum::PauliSum{TermType,CoeffType}, pstr::PauliString{TermType,CoeffType})

In-place subtraction a `PauliString` from a `PauliSum`. 
"""
function subtract!(psum::PauliSum{TermType,CoeffType}, pstr::PauliString{TermType,CoeffType}) where {TermType,CoeffType}
    _checknumberofqubits(psum, pstr)
    add!(psum.terms, pstr.term, -pstr.coeff)
    return psum
end

"""
    subtract!(psum1::PauliSum, psum2::PauliSum)

In-place subtraction a `PauliSum` from a `PauliSum`.
"""
function subtract!(psum1::PauliSum{TermType,CoeffType}, psum2::PauliSum{TermType,CoeffType}) where {TermType,CoeffType}
    _checknumberofqubits(psum1, psum2)
    subtract!(psum1.terms, psum2.terms)
    return psum1
end

"""
    subtract!(psum1::Dict{TermType,CoeffType}, psum2::Dict{TermType,CoeffType})

In-place subtraction a Pauli sum dict a Pauli sum dict.
"""
function subtract!(psum1::Dict{TermType,CoeffType}, psum2::Dict{TermType,CoeffType}) where {TermType,CoeffType}
    for (pstr, coeff) in psum2
        add!(psum1, pstr, -coeff)
    end
    return psum1
end

"""
    subtract!(psum::PauliSum{TermType, CoeffType}, pstr::TermType, coeff::CoeffType)

In-place subtraction of a Pauli string `pstr` with coefficient `coeff` to a `PauliSum`.
"""
function subtract!(psum::PauliSum{TermType,CoeffType}, pstr::TermType, coeff::CoeffType) where {TermType,CoeffType}
    subtract!(psum.terms, pstr, coeff)
    return psum
end

"""
    subtract!(psum::Dict{TermType, CoeffType}, pstr::TermType, coeff::CoeffType)

In-place subtraction of a Pauli string `pstr` with coefficient `coeff` to a Pauli sum dictionary.
"""
function subtract!(psum::Dict{TermType,CoeffType}, pstr::TermType, coeff::CoeffType) where {TermType,CoeffType}
    # don't add if the coefficient is 0
    if getnumcoeff(coeff) == 0
        return psum
    end

    if haskey(psum, pstr)
        # if the new coefficient is exactly 0, delete the key
        new_coeff = psum[pstr] - coeff
        if getnumcoeff(new_coeff) == 0
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
The type of the Pauli string needs to be the keytype of the dictionary, and the coefficient `coeff` needs to be the keytype.
"""
function set!(psum::PauliSum{TermType,CoeffType}, pstr::TermType, coeff::CoeffType) where {TermType,CoeffType}
    # TODO:Allow for truncation at this level
    set!(psum.terms, pstr, coeff)
    return psum
end

"""
    set!(psum::Dict{TermType, CoeffType}, pstr::TermType, coeff::CoeffType)

In-place setting the coefficient of a Pauli string in a Pauli sum dictionary.
The type of the Pauli string needs to be the keytype of the dictionary, and the coefficient `coeff` needs to be the keytype.
"""
function set!(psum::Dict{TermType,CoeffType}, pstr::TermType, coeff::CoeffType) where {TermType,CoeffType}
    # delete if the coefficient would be set to 0
    if getnumcoeff(coeff) == 0
        delete!(psum, pstr)

    else
        psum[pstr] = coeff
    end
    return psum
end

"""
    empty!(psum::PauliSum)

Empty the `PauliSum` by emtpying the dictionary on the `terms` fields. 
"""
Base.empty!(psum::PauliSum) = empty!(psum.terms)

## Helper functions
function similar(psum::PauliSum{TermType,CoeffType}) where {TermType,CoeffType}
    return PauliSum(psum.nqubits, similar(psum.terms))
end

function similar(psum::Dict{TermType,CoeffType}) where {TermType,CoeffType}
    new_psum = Dict{TermType,CoeffType}()
    sizehint!(new_psum, length(psum))
    return new_psum
end

"""
Converts coefficient to a float-type it it is an integer-type because the Pauli dictionaries need to be strictly typed and will likely become floats during propagation through a circuit.
"""
function _convertcoefficients(coeff)
    if isa(coeff, Integer)
        return Float64(coeff)
    elseif isa(coeff, Complex)
        return Complex{Float64}(coeff)
    else
        return coeff
    end
end

"""
Checks whether the number of qubits `nqubits` is the same between our datatypes.
"""
function _checknumberofqubits(nqubits::Int, pobj::Union{PauliString,PauliSum})
    if nqubits != pobj.nqubits
        throw(
            ArgumentError(
                "Number of qubits ($(nqubits)) must equal number of qubits ($(pobj.nqubits)) in $(typeof(pobj))"
            )
        )
    end
end

"""
Checks whether the number of qubits `nqubits` is the same between as the length of the vector `pstr`.
"""
function _checknumberofqubits(nqubits::Int, pstr)
    if nqubits != length(pstr)
        throw(
            ArgumentError(
                "Number of qubits ($(op1.nqubits)) must equal number of qubits ($(length(pstr))) in $(typeof(pstr))"
            )
        )
    end
end

"""
Checks whether the number of qubits `nqubits` is the same between our datatypes.
"""
function _checknumberofqubits(pobj1::Union{PauliString,PauliSum}, pobj2::Union{PauliString,PauliSum})
    if pobj1.nqubits != pobj2.nqubits
        throw(
            ArgumentError(
                "Number of qubits ($(pobj1.nqubits)) in $(typeof(pobj1)) must equal number of qubits ($(pobj2.nqubits)) in $(typeof(pobj2))"
            )
        )
    end
end