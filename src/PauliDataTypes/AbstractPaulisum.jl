"""
    AbstractPauliSum

Abstract type for objects represented sums of Paulis with coefficients.
"""
abstract type AbstractPauliSum end

nqubits(psum::AbstractPauliSum) = throw(ErrorException("nqubits not implemented for type $(typeof(psum))."))

"""
    iterate(psum::AbstractPauliSum)

Iterator interface for `AbstractPauliSum` types.
Yields tuples of `(paulis, coefficient)`.
"""
function Base.iterate(psum::AbstractPauliSum)
    # 1. Create the iterator we are delegating to
    iter = zip(paulis(psum), coefficients(psum))

    # 2. Start its iteration
    next = iterate(iter)

    # 3. Return the first item and a new state tuple: (iterator, iterator_state)
    #    We use a ternary operator for compactness.
    return next === nothing ? nothing : (next[1], (iter, next[2]))
end

function Base.iterate(psum::AbstractPauliSum, state)
    # 1. Unpack the state tuple
    (iter, inner_state) = state

    # 2. Continue the delegated iteration
    next = iterate(iter, inner_state)

    # 3. Return the next item and the updated state tuple
    return next === nothing ? nothing : (next[1], (iter, next[2]))
end


Base.length(psum::AbstractPauliSum) = length(paulis(psum))


function numcoefftype(psum::AbstractPauliSum)
    return numcoefftype(coefftype(psum))
end

function numcoefftype(CT::Type{<:Number})
    return CT
end


"""
    norm(psum::AbstractPauliSum, L=2)

Calculate the norm of a Pauli sum `psum` with respect to the `L`-norm. 
Calls `LinearAlgebra.norm(coefficients(psum))`.
"""
function LinearAlgebra.norm(psum::AbstractPauliSum, L::Real=2)
    if length(psum) == 0
        return zero(numcoefftype(psum))
    end
    return LinearAlgebra.norm((tonumber(coeff) for coeff in coefficients(psum)), L)
end