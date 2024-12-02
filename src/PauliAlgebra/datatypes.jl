import Base: *
import Base: /
import Base: +
import Base: -

"""
    PauliString(nqubits::Int, operator::OpType, coeff::CoeffType)

`PauliString` is a struct that represents a Pauli operator acting on `nqubits` qubits.
"""
struct PauliString{OpType<:PauliStringType,CoeffType}
    nqubits::Int
    operator::OpType  # TODO: rename
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

Constructor for a `PauliString` on `nqubits` qubits from a list of Symbols representing non-identity Pauli operator on qubits `qinds` with coefficient `coeff`.
"""
function PauliString(nqubits, pstr, qinds, coeff=1.0)
    pauli = symboltoint(nqubits, pstr, qinds)
    coeff = _convertcoefficients(coeff)
    return PauliString(nqubits, pauli, coeff)
end

"""
    paulitype(pstr::PauliString)

Get the Pauli integer type of a `PauliString`.
"""
function paulitype(pstr::PauliString)
    return typeof(pstr.operator)
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
    pauli_string = inttostring(pstr.operator, pstr.nqubits)
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
    PauliSum(nqubits::Int, op_dict::Dict)

`PauliSum`` is a struct that represents a sum of Pauli operators acting on `nqubits` qubits.
It is a wrapper around a dictionary Dict(Pauli string => coefficient}, where the Pauli strings are typically unsigned Integers for efficiency reasons.
"""
struct PauliSum{OpType<:Integer,CoeffType}
    nqubits::Int
    op_dict::Dict{OpType,CoeffType}
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
    OpType = getinttype(nq)
    return PauliSum(nq, Dict{OpType,ELTYPE}())
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
    PauliSum(nq::Integer, pstr::PauliString{OpType,CoeffType}) where {OpType,CoeffType}

Constructor for a `PauliSum` on `nqubits` qubits from a `PauliString`.
"""
function PauliSum(nq::Integer, pstr::PauliString{OpType,CoeffType}) where {OpType,CoeffType}
    _checknumberofqubits(nq, pstr)
    return PauliSum(nq, Dict{OpType,CoeffType}(pstr.operator => pstr.coeff))
end

"""
    paulis(psum::PauliSum)

Returns an iterator over the integer pauli strings of a `PauliSum`.
Call `topaulistrings` to receive entries as `PauliString`s.
"""
function paulis(psum::PauliSum)
    return keys(psum.op_dict)
end

"""
    coefficients(psum::PauliSum)

Returns an iterator over the coefficients of a `PauliSum`.
Call `topaulistrings` to receive entries as `PauliString`s.
"""
function coefficients(psum::PauliSum)
    return values(psum.op_dict)
end

"""
    paulitype(psum::PauliSum)

Get the Pauli integer type of a `PauliSum`.
"""
function paulitype(psum::PauliSum)
    return keytype(psum.op_dict)
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
    return valtype(psum.op_dict)
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
    return numcoefftype(psum.op_dict)
end

"""
    numcoefftype(psum::Dict)

Get the type of the numerical coefficient of a pauli sum dict. 
Will get the type of the `coeff` field of a potential PathProperties type.
"""
function numcoefftype(psum::Dict)
    return numcoefftype(first(values(psum)))
end

"""
    numcoefftype(::Number)

Get the type of the number.
"""
function numcoefftype(::T) where {T<:Number}
    return T
end

"""
    getcoeff(psum::PauliSum{OpType,CoeffType}, operator::Integer) where {OpType,CoeffType}

Get the coefficient of an integer Pauli string in a `PauliSum`. Defaults to 0 if the Pauli string is not in the `PauliSum`.
"""
function getcoeff(psum::PauliSum{OpType,CoeffType}, operator::Integer) where {OpType,CoeffType}
    return get(psum.op_dict, operator, CoeffType(0))
end

"""
    getcoeff(psum::PauliSum{OpType,CoeffType1}, pstr::PauliString{OpType,CoeffType2}) where {OpType,CoeffType1,CoeffType2}

Get the coefficient of a `PauliString` in a `PauliSum`. Defaults to 0 if the Pauli string is not in the `PauliSum`.
"""
function getcoeff(psum::PauliSum{OpType,CoeffType1}, pstr::PauliString{OpType,CoeffType2}) where {OpType,CoeffType1,CoeffType2}
    return get(psum.op_dict, pstr.operator, CoeffType1(0))
end


"""
    getcoeff(psum::PauliSum{OpType,CoeffType}, pauli::Symbol, qind::Integer) where {OpType,CoeffType}

Get the coefficient of a Pauli string in a `PauliSum` by providing the Pauli string as a Symbol acting on qubit `qind`. 
This is consistent with how operators can be added to a `PauliSum` via `add!()`. Defaults to 0 if the operator is not in the `PauliSum`.
"""
function getcoeff(psum::PauliSum{OpType,CoeffType}, pauli::Symbol, qind::Integer) where {OpType,CoeffType}
    return getcoeff(psum, symboltoint(psum.nqubits, pauli, qind))
end

"""
    getcoeff(psum::PauliSum{OpType,CoeffType}, pstr, qinds) where {OpType,CoeffType}

Get the coefficient of a Pauli string in a `PauliSum` by providing the Pauli string `pstr` as a vector of Symbols acting on qubits `qinds`. 
This is consistent with how operators can be added to a `PauliSum` via `add!()`. Defaults to 0 if the operator is not in the `PauliSum`.
"""
function getcoeff(psum::PauliSum{OpType,CoeffType}, pstr, qinds) where {OpType,CoeffType}
    return getcoeff(psum, symboltoint(psum.nqubits, pstr, qinds))
end

"""
    getcoeff(psum::PauliSum{OpType,CoeffType}, pstr) where {OpType,CoeffType}

Get the coefficient of a Pauli string in a `PauliSum` by providing the Pauli string `pstr` as a vector of Symbols acting on all qubits. 
This is consistent with how operators can be added to a `PauliSum` via `add!()`. Defaults to 0 if the operator is not in the `PauliSum`.
"""
function getcoeff(psum::PauliSum{OpType,CoeffType}, pstr) where {OpType,CoeffType}
    return getcoeff(psum, symboltoint(pstr))
end

"""
    getnumcoeff(val::Number)

Trivial function returning a numerical value of a number.
"""
function getnumcoeff(val::Number)
    return val
end

# TODO: Add functions for extracting paulis and coefficients from the PauliSum (as iterable)

"""
    topaulistrings(psum::PauliSum)

Returns the Pauli strings in a `PauliSum` and their coefficients as a list of `PauliString`.
"""
topaulistrings(psum::PauliSum) = [PauliString(psum.nqubits, pauli, coeff) for (pauli, coeff) in psum.op_dict]

"""
Copy a `PauliSum` by copying its `op_dict`.
"""
Base.copy(psum::PauliSum) = PauliSum(psum.nqubits, copy(psum.op_dict))

"""
Iterator for `PauliSum` returns the iterator over its `op_dict`.
"""
Base.iterate(psum::PauliSum, state=1) = iterate(psum.op_dict, state)


import Base.show
"""
Pretty print for `PauliSum`.
"""
function show(io::IO, psum::PauliSum)
    if length(psum.op_dict) == 0
        dict_string = "(no operators)"
    else
        dict_string = _getprettystr(psum.op_dict, psum.nqubits)
    end
    print(io, "PauliSum(nqubits: $(psum.nqubits), $dict_string)")
end

import Base.length
"""
    length(psum::PauliSum)

Number of terms in the `PauliSum`.
"""
length(psum::PauliSum) = length(psum.op_dict)

import Base: ==
"""
    ==(psum1::PauliSum, psum2::PauliSum)

Equality check for `PauliSum`.
"""
function ==(psum1::PauliSum, psum2::PauliSum)
    if psum1.nqubits != psum2.nqubits
        return false
    end

    return psum1.op_dict == psum2.op_dict
end

"""
    mult!(psum::PauliSum, c::Number)

Multiply a `PauliSum` by a scalar `c` in-place.
"""
function mult!(psum::PauliSum, c::Number)
    # multiply in-place
    for (k, v) in psum.op_dict
        psum.op_dict[k] *= c
    end
    return psum
end

"""
    *(psum::PauliSum, c::Number)

Multiply a `PauliSum` by a scalar `c`. This copies the PauliSum.
"""
function *(psum::PauliSum, c::Number)
    ps_copy = copy(psum)  #TODO: make sure deepcopy is not needed
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
    ps_copy = copy(psum)  #TODO: make sure deepcopy is not needed
    mult!(ps_copy, 1 / c)
    return ps_copy
end

"""
    +(pstr1::PauliString, pstr2::PauliString)

Addition of two `PauliString`s. Returns a PauliSum.
"""
function +(pstr1::PauliString, pstr2::PauliString)
    psum = PauliSum(pstr1) # or deepcopy?
    add!(psum, pstr2)
    return psum
end

"""
    +(psum::PauliSum, pstr::PauliString)

Addition of a `PauliString` to a `PauliSum`. Returns a `PauliSum`.
"""
function +(psum::PauliSum, pstr::PauliString)
    psum = copy(psum) # or deepcopy?
    add!(psum, pstr)
    return psum
end

"""
    +(psum1::PauliSum, psum2::PauliSum)
Addition of two `PauliSum`s. Returns a `PauliSum`.
"""
function +(psum1::PauliSum, psum2::PauliSum)
    psum1 = copy(psum1) # or deepcopy?
    add!(psum1, psum2)
    return psum1
end

"""
    add!(psum::PauliSum, pstr::PauliString)

Addition of a `PauliString` to a `PauliSum`. Changes the `PauliSum` in-place.
"""
function add!(psum::PauliSum, pstr::PauliString)
    _checknumberofqubits(psum, pstr)
    psum.op_dict[pstr.operator] = get(psum.op_dict, pstr.operator, keytype(psum.op_dict)(0.0)) + pstr.coeff
    return psum
end

"""
    add!(psum1::PauliSum, psum2::PauliSum; precision=_DEFAULT_PRECISION)

Addition of two `PauliSum`s. Changes the first `PauliSum` in-place.
Uses a default precision for coefficients under which a coefficient is considered to be 0.
"""
function add!(psum1::PauliSum, psum2::PauliSum, precision=_DEFAULT_PRECISION)
    _checknumberofqubits(psum1, psum2)
    for (operator, coeff) in psum2.op_dict
        if haskey(psum1.op_dict, operator)
            psum1.op_dict[operator] += coeff

            # Remove the operator if the resulting coefficient is small
            if abs(psum1.op_dict[operator]) < precision
                delete!(psum1.op_dict, operator)
            end

        else
            psum1.op_dict[operator] = -coeff
        end
    end

    return psum1
end

"""
    add!(psum::PauliSum, pauli::Symbol, qind::Integer, coeff=1.0)

In-place addition a Pauli string to `PauliSum` by providing the Pauli string as a Symbol acting on qubit `qind`.
Coefficient defaults to 1.0.
"""
function add!(psum::PauliSum, pauli::Symbol, qind::Integer, coeff=1.0)
    return add!(psum, PauliString(psum.nqubits, pauli, qind, coeff))
end


# TODO: add! for PauliStringType with coefficient. Use these functions throughout `propagate`

"""
    add!(psum::PauliSum, pstr, qinds, coeff=1.0)

In-place addition a Pauli string to `PauliSum` by providing the Pauli string `pstr` as a vector of Symbols acting on qubits `qinds`.
Coefficient defaults to 1.0.
"""
function add!(psum::PauliSum, pstr, qinds, coeff=1.0) # TODO: don't strictly type qinds or similar here or elsewhere
    return add!(psum, PauliString(psum.nqubits, pstr, qinds, coeff))
end

## Substraction
const _DEFAULT_PRECISION = 1e-12
"""
    -(pstr1::PauliString, pstr2::PauliString)

Subtraction of two `PauliString`s. Returns a PauliSum. Uses the default precision of subtract!().
"""
function -(pstr1::PauliString, pstr2::PauliString)
    psum = PauliSum(pstr1) # or deepcopy?
    subtract!(psum, pstr2)
    return psum
end

"""
    -(psum::PauliSum, pstr::PauliString)

Subtraction of a `PauliString` to a `PauliSum`. Returns a `PauliSum`. Uses the default precision of subtract!().
"""
function -(psum::PauliSum, pstr::PauliString)
    psum = copy(psum) # or deepcopy?
    subtract!(psum, pstr)
    return psum
end

"""
    -(psum1::PauliSum, psum2::PauliSum)

Subtract of two `PauliSum`s. Returns a `PauliSum`. Uses the default precision of subtract!().
"""
function -(psum1::PauliSum, psum2::PauliSum)
    psum1 = copy(psum1) # or deepcopy?
    subtract!(psum1, psum2)
    return psum1
end

"""
    ubtract!(psum::PauliSum, pstr::PauliString; precision=_DEFAULT_PRECISION)

In-place subtraction a `PauliString` from a `PauliSum`. 
Uses a default precision for coefficients under which a coefficient is considered to be 0.
"""
function subtract!(psum::PauliSum, pstr::PauliString; precision=_DEFAULT_PRECISION)
    _checknumberofqubits(psum, pstr)
    if haskey(psum.op_dict, pstr.operator)
        psum.op_dict[pstr.operator] -= pstr.coeff

        # Remove the operator if the resulting coefficient is small
        if abs(psum.op_dict[pstr.operator]) < precision
            delete!(psum.op_dict, pstr.operator)
        end

    else
        psum.op_dict[pstr.operator] = -pstr.coeff
    end
    return psum
end

"""
    subtract!(psum1::PauliSum, psum2::PauliSum; precision=_DEFAULT_PRECISION)

In-place subtraction a `PauliSum` from a `PauliSum`. Uses a default precision for coefficients under which a coefficient is considered to be 0.
"""
function subtract!(psum1::PauliSum, psum2::PauliSum; precision=_DEFAULT_PRECISION)
    _checknumberofqubits(psum1, psum2)
    for (operator, coeff) in psum2.op_dict
        if haskey(psum1.op_dict, operator)
            psum1.op_dict[operator] -= coeff

            # Remove the operator if the resulting coefficient is small
            if abs(psum1.op_dict[operator]) < precision
                delete!(psum1.op_dict, operator)
            end

        else
            psum1.op_dict[operator] = -coeff
        end
    end

    return psum1
end

## Helper functions
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