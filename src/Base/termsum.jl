
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
terms(term_sum::TS) where TS<:AbstractTermSum = _thrownotimplemented(TS, :terms)
coefficients(term_sum::TS) where TS<:AbstractTermSum = _thrownotimplemented(TS, :coefficients)
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
    if haskey(term_dict, trm)
        return term_dict[trm]
    else
        return zero(coefftype(term_sum))
    end
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


### Default functions defined for all TermSum types
function Base.iterate(term_sum::AbstractTermSum)
    # 1. Create the iterator we are delegating to
    iter = zip(terms(term_sum), coefficients(term_sum))

    # 2. Start its iteration
    next = iterate(iter)

    # 3. Return the first item and a new state tuple: (iterator, iterator_state)
    #    We use a ternary operator for compactness.
    return next === nothing ? nothing : (next[1], (iter, next[2]))
end

function Base.iterate(term_sum::AbstractTermSum, state)
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

# TODO: make general:
Base.empty!(term_sum::TS) where TS<:AbstractTermSum = empty!(storage(term_sum))
Base.delete!(term_sum::TS, term) where TS<:AbstractTermSum = delete!(storage(term_sum), term)

"""
    similar(term_sum::AbstractTermSum)
    
Create a similar, empty TermSum of the same type as `term_sum`.
By default, this is implemented via `deepcopy` and `empty!`.
Consider overloading for performance if needed.
"""
function Base.similar(term_sum::TS) where TS<:AbstractTermSum
    similar_term_sum = deepcopy(term_sum)
    empty!(similar_term_sum)
    return similar_term_sum
end


### Adding, setting, and deleting
"""
    add!(term_sum::AbstractTermSum, term, coeff)

Add `coeff` to the coefficient of `term` in `term_sum`.
Calls `add!(storage(term_sum), term, coeff)` internally.
For custom behavior, overload `storage()` and/or `add!` for the specific TermSum type.
"""
function add!(term_sum::TS, term, coeff) where TS<:AbstractTermSum
    add!(storage(term_sum), term, coeff)
    return term_sum
end


function add!(dict_storage::Dict, term, coeff)
    if haskey(dict_storage, term)
        dict_storage[term] += coeff
    else
        dict_storage[term] = coeff
    end
    return dict_storage
end

function add!(vector_storage::Tuple{<:AbstractVector,<:AbstractVector}, term, coeff)
    terms_vec, coeffs_vec = vector_storage
    for i in eachindex(terms_vec)
        if terms_vec[i] == term
            coeffs_vec[i] += coeff
            return vector_storage
        end
    end
    push!(terms_vec, term)
    push!(coeffs_vec, coeff)
    return vector_storage
end

# """
#     add!(term_sum1::AbstractTermSum, term_sum2::AbstractTermSum)

# Add all terms and coefficients from `term_sum2` into `term_sum1` via `mergewith!(+, ...)`.
# """
add!(term_sum1::AbstractTermSum, term_sum2::AbstractTermSum) = _thrownotimplemented(typeof(term_sum1), :add!)
# TODO: implement this via merge!()


"""
    set!(term_sum::AbstractTermSum, term, coeff)
    
Set the coefficient of `term` in `term_sum` to `coeff`.
Calls `set!(storage(term_sum), term, coeff)` internally.
For custom behavior, overload `storage()` and/or `set!` for the specific TermSum type.
"""
set!(term_sum::TS, term, coeff) where TS<:AbstractTermSum = set!(storage(term_sum), term, coeff)

"""
    set!(term_storage::Dict, term, coeff)

Default implementation of `set!()` for `Dict` storage.
"""
function set!(dict_storage::Dict, term, coeff)
    dict_storage[term] = coeff
    return dict_storage
end

function set!(vector_storage::Tuple{<:AbstractVector,<:AbstractVector}, term, coeff)
    terms_vec, coeffs_vec = vector_storage
    for i in eachindex(terms_vec)
        if terms_vec[i] == term
            coeffs_vec[i] = coeff
            return vector_storage
        end
    end
    push!(terms_vec, term)
    push!(coeffs_vec, coeff)
    return vector_storage
end

"""
    mult!(term_sum::AbstractTermSum, scalar::Number)

Multiply all coefficients in `term_sum` by `scalar`.
Calls `mult!(storage(term_sum), scalar)` internally.
For custom behavior, overload `storage()` and/or `mult!` for the specific TermSum type.
"""
function mult!(term_sum::TS, scalar::Number) where TS<:AbstractTermSum
    return mult!(storage(term_sum), scalar)
end


function mult!(dict_storage::Dict, scalar::Number)
    for (term, coeff) in dict_storage
        dict_storage[term] = coeff * scalar
    end
    return dict_storage
end

function mult!(vector_storage::Tuple{<:AbstractVector,<:AbstractVector}, scalar::Number)
    terms_vec, coeffs_vec = vector_storage
    coeffs_vec .*= scalar
    return vector_storage
end


### Short out-of-place algebra
# TODO all of these. Requires a proper add/merge implementation


function Base.:+(term_sum1::TS, term_sum2::TS) where TS<:AbstractTermSum
    result = deepcopy(term_sum1)
    mergewith!(+, result, term_sum2)
    return result
end

function Base.:-(term_sum1::TS, term_sum2::TS) where TS<:AbstractTermSum
    result = deepcopy(term_sum1)
    mergewith!(-, result, term_sum2)
    return result
end

function Base.:*(term_sum::TS, scalar) where TS<:AbstractTermSum
    result = deepcopy(term_sum)
    mult!(result, scalar)
    return result
end


# check for equality by equaility on all fields
function Base.:(==)(a::AbstractTermSum, b::AbstractTermSum)
    return all(field -> getfield(a, field) == getfield(b, field), fieldnames(typeof(a)))
end