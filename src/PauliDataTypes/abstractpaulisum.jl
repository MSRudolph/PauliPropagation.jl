"""
    AbstractPauliSum <: AbstractTermSum

Abstract type for objects represented sums of Paulis with coefficients.
"""
abstract type AbstractPauliSum <: AbstractTermSum end

nqubits(psum::AbstractPauliSum) = throw(ErrorException("nqubits() not implemented for type $(typeof(psum))."))
PropagationBase.nsites(psum::AbstractPauliSum) = nqubits(psum)

"""
    paulis(psum::AbstractPauliSum)

Returns an iterator over the integer pauli strings of an `AbstractPauliSum`.
Call `topaulistrings` to receive entries as `PauliString`s.
"""
paulis(psum::AbstractPauliSum) = terms(psum)

"""
    coefficients(psum::AbstractPauliSum)

Returns an iterator over the coefficients of a `PauliSum`.
Call `topaulistrings` to receive entries as `PauliString`s.
"""
PropagationBase.coefficients

"""
    paulitype(psum::AbstractPauliSum)

Get the Pauli integer type of a `AbstractPauliSum` object.
"""
paulitype(psum::AbstractPauliSum) = eltype(paulis(psum))


"""
    topaulistrings(psum::AbstractPauliSum)

Returns the Pauli strings in a, `AbstractPauliSum` and their coefficients as a list of `PauliString`.
"""
topaulistrings(psum::AbstractPauliSum) = [PauliString(psum.nqubits, pauli, coeff) for (pauli, coeff) in zip(paulis(psum), coefficients(psum))]


## The symbol conversions for getcoeff()

"""
    getcoeff(psum::AbstractPauliSum, pauli::Symbol, qind::Integer)

Get the coefficient of a Pauli string in an `AbstractPauliSum` by providing the Pauli string as a Symbol acting on qubit `qind`. 
This is consistent with how Pauli strings can be added to a `PauliSum` via `add!()`. 
Defaults to 0 if the Pauli string is not in the `AbstractPauliSum`.
"""
function PropagationBase.getcoeff(psum::AbstractPauliSum, pauli::Symbol, qind::Integer)
    return getcoeff(psum, symboltoint(psum.nqubits, pauli, qind))
end

"""
    getcoeff(psum::AbstractPauliSum, pstr::Vector{Symbol}, qinds::Vector{Int})

Get the coefficient of a Pauli string in an `AbstractPauliSum` by providing the Pauli string `pstr` as a vector of Symbols acting on qubits `qinds`. 
This is consistent with how Pauli strings can be added to a `PauliSum` via `add!()`. 
Defaults to 0 if the Pauli string is not in the `AbstractPauliSum`.
"""
function PropagationBase.getcoeff(psum::AbstractPauliSum, pstr, qinds)
    return getcoeff(psum, symboltoint(psum.nqubits, pstr, qinds))
end

"""
    getcoeff(psum::AbstractPauliSum, pstr::Vector{Symbol})

Get the coefficient of a Pauli string in a `AbstractPauliSum` by providing the Pauli string `pstr` as a vector of Symbols acting on all qubits. 
This is consistent with how Pauli strings can be added to a `PauliSum` via `add!()`. 
Defaults to 0 if the Pauli string is not in the `AbstractPauliSum`.
"""
function PropagationBase.getcoeff(psum::AbstractPauliSum, pstr::Vector{Symbol})
    return getcoeff(psum, symboltoint(pstr))
end


"""
    filter!(filterfunc::Function, psum::AbstractPauliSum)

Filter a `AbstractPauliSum` by copying and removing all Pauli strings for which `filterfunc(pstr, coeff)` returns `false`.
"""
Base.filter(filterfunc::F, psum::AbstractPauliSum) where {F<:Function} = truncate!((pstr, coeff) -> !filterfunc(pstr, coeff), deepcopy(psum))

"""
    filter!(filterfunc::Function, psum::AbstractPauliSum)

Filter a `AbstractPauliSum` in-place by removing all Pauli strings for which `filterfunc(pstr, coeff)` returns `false`.
"""
filter!(filterfunc::F, psum::AbstractPauliSum) where {F<:Function} = truncate!((pstr, coeff) -> !filterfunc(pstr, coeff), psum)