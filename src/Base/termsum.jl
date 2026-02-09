
# TermSum is an abstract type for container types carrying terms and coefficients
# we expect that it 
abstract type AbstractTermSum end

### Interface functions to be defined for normal use of TermSum types 


abstract type StorageType end
struct DictStorage <: StorageType end
struct ArrayStorage <: StorageType end

function StorageType(term_sum::AbstractTermSum)
    if storage(term_sum) isa Dict
        return DictStorage()
    elseif storage(term_sum) isa Tuple{<:AbstractArray,<:AbstractArray}
        return ArrayStorage()
    else
        return _thrownotimplemented(typeof(term_sum), :StorageType)
    end
end

# storage() is expected to return the internal storage representation of the TermSum
# often this is a Dict{TermType,CoeffType} but it can be anything
storage(term_sum::TS) where TS<:AbstractTermSum = _thrownotimplemented(TS, :storage)

Base.length(term_sum::AbstractTermSum) = length(terms(term_sum))

nsites(term_sum::TS) where TS<:AbstractTermSum = _thrownotimplemented(TS, :nsites)

terms(::StorageType, term_sum::TS) where TS<:AbstractTermSum = _thrownotimplemented(TS, :terms)
terms(term_sum::AbstractTermSum) = terms(StorageType(term_sum), term_sum)
terms(::DictStorage, term_sum::AbstractTermSum) = keys(storage(term_sum))
terms(::ArrayStorage, term_sum::AbstractTermSum) = storage(term_sum)[1]


coefficients(::StorageType, term_sum::TS) where TS<:AbstractTermSum = _thrownotimplemented(TS, :coefficients)
coefficients(term_sum::AbstractTermSum) = coefficients(StorageType(term_sum), term_sum)
coefficients(::DictStorage, term_sum::AbstractTermSum) = values(storage(term_sum))
coefficients(::ArrayStorage, term_sum::AbstractTermSum) = storage(term_sum)[2]
coeffs(term_sum::AbstractTermSum) = coefficients(term_sum)

# receives the object
termtype(term_sum::TS) where TS<:AbstractTermSum = eltype(terms(term_sum))
coefftype(term_sum::TS) where TS<:AbstractTermSum = eltype(coefficients(term_sum))

# this is used to determine type-stable return values for numerical operations
numcoefftype(term_sum::TS) where TS<:AbstractTermSum = numcoefftype(coefftype(term_sum))
numcoefftype(::Type{CT}) where CT<:Number = CT
numcoefftype(::Type{T}) where T = _thrownotimplemented(T, :numcoefftype)


function getcoeff(term_sum::AbstractTermSum, trm)
    getcoeff(StorageType(term_sum), term_sum, trm)
end

function getcoeff(::Type{DictStorage}, term_sum::AbstractTermSum, trm)
    term_dict = storage(term_sum)
    return get(term_dict, trm, zero(coefftype(term_sum)))
end

# default implementation
function getcoeff(::ST, term_sum::AbstractTermSum, trm) where {ST<:StorageType}
    # TODO: GPU kernel for this
    val = zero(coefftype(term_sum))
    for (term, coeff) in zip(terms(term_sum), coefficients(term_sum))
        if term == trm
            val += coeff
        end
    end
    return val
end

# this assumes everything is merged and de-duplicated
# may result in wrong results if not
function getmergedcoeff(term_sum::AbstractTermSum, trm)
    return getmergedcoeff(StorageType(term_sum), term_sum, trm)
end

# for the DictStorage, this is the same as getcoeff
function getmergedcoeff(::DictStorage, term_sum::AbstractTermSum, trm)
    return getcoeff(DictStorage(), term_sum, trm)
end

# everywhere else, we do a linear search
function getmergedcoeff(::StorageType, term_sum::AbstractTermSum, trm)
    terms_vec, coeffs_vec = storage(term_sum)
    for (term, coeff) in zip(terms_vec, coeffs_vec)
        if term == trm
            return coeff
        end
    end
    return zero(coefftype(term_sum))
end


### Default functions defined for all TermSum types
@inline function Base.iterate(term_sum::AbstractTermSum, args...)
    return _iterate(StorageType(term_sum), term_sum, args...)
end

@inline function _iterate(::DictStorage, term_sum::AbstractTermSum, args...)
    dict_storage = storage(term_sum)
    return iterate(dict_storage, args...)
end


@inline function _iterate(::StorageType, term_sum::AbstractTermSum)
    # 1. Create the iterator we are delegating to
    iter = zip(terms(term_sum), coefficients(term_sum))

    # 2. Start its iteration
    next = iterate(iter)

    # 3. Return the first item and a new state tuple: (iterator, iterator_state)
    #    We use a ternary operator for compactness.
    return next === nothing ? nothing : (next[1], (iter, next[2]))
end

@inline function _iterate(::StorageType, term_sum::AbstractTermSum, state)
    # 1. Unpack the state tuple
    (iter, inner_state) = state

    # 2. Continue the delegated iteration
    next = iterate(iter, inner_state)

    # 3. Return the next item and the updated state tuple
    return next === nothing ? nothing : (next[1], (iter, next[2]))
end


"""
    norm(psum::AbstractTermSum, L=2)

Calculate the norm of a Pauli sum `psum` with respect to the `L`-norm. 
Calls `LinearAlgebra.norm(coefficients(psum))`.
"""
function LinearAlgebra.norm(psum::AbstractTermSum, L::Real=2)
    if length(psum) == 0
        return zero(numcoefftype(psum))
    end
    return LinearAlgebra.norm((tonumber(coeff) for coeff in coefficients(psum)), L)
end


### Adding, setting, and deleting
"""
    add!(term_sum::AbstractTermSum, term, coeff)

Add `coeff` to the coefficient of `term` in `term_sum`.
Calls `_add!(StorageType(term_sum), term_sum, term, coeff)` internally.
For custom behavior, overload `storage()` and/or `_add!` for the specific TermSum type.
"""
@inline function add!(term_sum::AbstractTermSum, term, coeff)
    _add!(StorageType(term_sum), term_sum, term, coeff)
    return term_sum
end


@inline function _add!(::DictStorage, term_sum::AbstractTermSum, term, coeff)
    dict_storage = storage(term_sum)
    if haskey(dict_storage, term)
        dict_storage[term] += coeff
    else
        dict_storage[term] = coeff
    end
    return term_sum
end

@inline function _add!(::ArrayStorage, term_sum::AbstractTermSum, term, coeff)
    terms_vec, coeffs_vec = storage(term_sum)
    ind = findfirst(t -> t == term, terms_vec)
    if !isnothing(ind)
        coeffs_vec[ind] += coeff
    else
        push!(terms_vec, term)
        push!(coeffs_vec, coeff)
    end
    return term_sum
end

@inline function _add!(::StorageType, term_sum::AbstractTermSum, term, coeff)
    thrownotimplemented(typeof(term_sum), :_add!)
end


function add!(term_sum1::AbstractTermSum, term_sum2::AbstractTermSum)
    for (term, coeff) in term_sum2
        add!(term_sum1, term, coeff)
    end
    return term_sum1
end

function subtract!(term_sum1::AbstractTermSum, term_sum2::AbstractTermSum)
    for (term, coeff) in term_sum2
        add!(term_sum1, term, -coeff)
    end
    return term_sum1
end


"""
    set!(term_sum::AbstractTermSum, term, coeff)
    
Set the coefficient of `term` in `term_sum` to `coeff`.
Calls `_set!(StorageType(term_sum), term_sum, term, coeff)` internally.
For custom behavior, overload `storage()` and/or `_set!` for the specific TermSum type.
"""
@inline function set!(term_sum::AbstractTermSum, term, coeff)
    _set!(StorageType(term_sum), term_sum, term, coeff)
    return term_sum
end


@inline function _set!(::DictStorage, term_sum::AbstractTermSum, term, coeff)
    dict_storage = storage(term_sum)
    dict_storage[term] = coeff
    return term_sum
end

@inline function _set!(::ArrayStorage, term_sum::AbstractTermSum, term, coeff)
    terms_vec, coeffs_vec = storage(term_sum)
    ind = findfirst(t -> t == term, terms_vec)
    if !isnothing(ind)
        coeffs_vec[ind] = coeff
    else
        push!(terms_vec, term)
        push!(coeffs_vec, coeff)
    end
    return term_sum
end

@inline function _set!(ST::StorageType, term_sum::AbstractTermSum, term, coeff)
    old_coeff = getcoeff(term_sum, term)
    if old_coeff == zero(coefftype(term_sum))
        _add!(ST, term_sum, term, coeff)
    else
        delta = coeff - old_coeff
        _add!(ST, term_sum, term, delta)
    end
    return term_sum
end

"""
    mult!(term_sum::AbstractTermSum, scalar::Number)

Multiply all coefficients in `term_sum` by `scalar`.
Calls `mult!(StorageType(term_sum), term_sum, scalar)` internally.
For custom behavior, overload `storage()` and/or `mult!` for the specific TermSum type.
"""
function mult!(term_sum::AbstractTermSum, scalar::Number)
    return mult!(StorageType(term_sum), term_sum, scalar)
end


function mult!(::DictStorage, term_sum::AbstractTermSum, scalar::Number)
    dict_storage = storage(term_sum)
    for (term, coeff) in dict_storage
        dict_storage[term] = coeff * scalar
    end
    return term_sum
end

function mult!(::ArrayStorage, term_sum::AbstractTermSum, scalar::Number)
    terms_vec, coeffs_vec = storage(term_sum)
    coeffs_vec .*= scalar
    return term_sum
end

# super slow default
function mult!(::StorageType, term_sum::AbstractTermSum, scalar::Number)
    for (term, coeff) in zip(terms(term_sum), coefficients(term_sum))
        set!(term_sum, term, coeff * scalar)
    end
    return term_sum
end

function Base.delete!(term_sum::AbstractTermSum, term)
    delete!(StorageType(term_sum), term_sum, term)
    return term_sum
end

function Base.delete!(::DictStorage, term_sum::AbstractTermSum, term)
    dict_storage = storage(term_sum)
    delete!(dict_storage, term)
    return term_sum
end

function Base.delete!(::ArrayStorage, term_sum::AbstractTermSum, term)
    terms_vec, coeffs_vec = storage(term_sum)
    ind = findfirst(t -> t == term, terms_vec)
    if !isnothing(ind)
        deleteat!(terms_vec, ind)
        deleteat!(coeffs_vec, ind)
    end
    return term_sum
end

function Base.delete!(::StorageType, term_sum::AbstractTermSum, term)
    # by default we set the coefficient to zero
    set!(term_sum, term, zero(coefftype(term_sum)))
    return term_sum
end

function Base.empty!(term_sum::AbstractTermSum)
    return empty!(StorageType(term_sum), term_sum)
end

function Base.empty!(::DictStorage, term_sum::AbstractTermSum)
    dict_storage = storage(term_sum)
    empty!(dict_storage)
    return term_sum
end

function Base.empty!(::ArrayStorage, term_sum::AbstractTermSum)
    terms_vec, coeffs_vec = storage(term_sum)
    empty!(terms_vec)
    empty!(coeffs_vec)
    return term_sum
end

function Base.empty!(::StorageType, term_sum::AbstractTermSum)
    for term in terms(term_sum)
        delete!(term_sum, term)
    end
    return term_sum
end


function Base.similar(term_sum::AbstractTermSum)
    similar_term_sum = deepcopy(term_sum)
    empty!(similar_term_sum)
    return similar_term_sum
end


### Short out-of-place algebra

function Base.:+(term_sum1::AbstractTermSum, term_sum2::AbstractTermSum)
    result = deepcopy(term_sum1)
    add!(result, term_sum2)
    return result
end


function Base.:-(term_sum1::AbstractTermSum, term_sum2::AbstractTermSum)
    result = deepcopy(term_sum1)
    subtract!(result, term_sum2)
    return result
end

function Base.:*(term_sum::AbstractTermSum, scalar)
    result = deepcopy(term_sum)
    mult!(result, scalar)
    return result
end

Base.:*(c::Number, term_sum::AbstractTermSum) = term_sum * c

Base.:/(term_sum::AbstractTermSum, scalar) = term_sum * (one(scalar) / scalar)


# check for equality by equaility on all fields
function Base.isequal(term_sum1::AbstractTermSum, term_sum2::AbstractTermSum)
    typeof(term_sum1) == typeof(term_sum2) || return false
    termtype(term_sum1) == termtype(term_sum2) || return false
    coefftype(term_sum1) == coefftype(term_sum2) || return false
    length(term_sum1) == length(term_sum2) || return false

    # call isequal on all fields of the types 
    return all(isequal(getfield(term_sum1, f), getfield(term_sum2, f)) for f in fieldnames(typeof(term_sum1)))
end

Base.:(==)(term_sum1::AbstractTermSum, term_sum2::AbstractTermSum) = isequal(term_sum1, term_sum2)

function Base.isapprox(term_sum1::AbstractTermSum, term_sum2::AbstractTermSum; approx_kwargs...)
    termtype(term_sum1) == termtype(term_sum2) || return false
    coefftype(term_sum1) == coefftype(term_sum2) || return false
    length(term_sum1) == length(term_sum2) || return false

    # call isapprox on all fields of the types
    return all(isapprox(getfield(term_sum1, f), getfield(term_sum2, f); approx_kwargs...) for f in fieldnames(typeof(term_sum1)))
end

# Base.:≈(term_sum1::AbstractTermSum, term_sum2::AbstractTermSum) = isapprox(term_sum1, term_sum2)