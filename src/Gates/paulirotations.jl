### paulirotations.jl
##
# A file for Pauli rotation gates, i.e., gates that are generated by a Pauli string.
# For example PauliRotation(:X, 2)  ->  RX_2(θ) = exp(-i θ/2 X_2) is a Pauli rotation gate.
##
###

# TODO: We should make these type-stable by using tuples and not vectors
"""
A type for a Pauli rotation gate carrying the gate generator and the qubit indices on which it acts.
"""
struct PauliRotation <: ParametrizedGate
    symbols::Vector{Symbol}
    qinds::Vector{Int}

    @doc """
        PauliRotation(symbols::Vector{Symbol}, qinds::Vector{Int})
        PauliRotation(symbol::Symbol, qind::Int)

    A parametrized Pauli rotation generated by the Pauli string `symbols` acting on the qubits `qinds`.
    For example PauliRotation(:X, 2) or PauliRotation([:X, :Y], [1, 2]).
    We follow the convention that those rotation are equivalent by conjugation of 
    ```math
    e^{-i θ/2 P}
    ```
    where `P` is the corresponding Pauli operator.
    """
    function PauliRotation(symbols, qinds)

        # turn symbols into vectors
        if isa(symbols, Symbol)
            symbols = [symbols]
        else
            symbols = vec(collect(symbols))
        end

        # check that the symbols are valid Pauli symbols
        if any(s -> s ∉ (:I, :X, :Y, :Z), symbols)
            throw(ArgumentError("Symbols must be `:I`, `:X`, `:Y`, or `:Z`. Got $symbols."))
        end

        # turn qinds into vectors
        qinds = vec(collect(qinds))
        # check that the qubit indices are positive integers
        _qinds_check(qinds)

        # check that the number of symbols matches the number of qubits
        if length(symbols) != length(qinds)
            throw(ArgumentError(
                "The number of symbols must match the number of qubits. " *
                "Got $(length(symbols)) symbols and $(length(qinds)) qubits."
            ))
        end

        return new(symbols, qinds)
    end
end


"""
    PauliRotation(symbols, qinds, theta)

Constructor for a frozen `PauliRotation` generated by the Pauli string `symbols` acting on the qubits `qinds`,
and with fixed parameter `theta`.
"""
function PauliRotation(symbols, qinds, theta)
    return FrozenGate(PauliRotation(symbols, qinds), theta)
end


"""
    tomatrix(gate::PauliRotation, theta)

Compute the unitary matrix for the `PauliRotation` gate with parameter `theta` in the computational 0/1 basis.
This is done by computing the matrix `U = cos(θ/2) I - i sin(θ/2) P` where `P` is the Pauli matrix corresponding to the `symbols`.
The returned unitary is returned in Schrödinger picture form. 
"""
function tomatrix(gate::PauliRotation, theta)

    nqubits = length(gate.qinds)
    pauli_basis_vec = getpaulibasis(nqubits)

    # The indices of the pauli matrix are sorted in ascending order
    sorted_indices = sortperm(gate.qinds)
    pauli = gate.symbols[sorted_indices]

    # These pauli matrices are normalized by sqrt(2)^nqubits
    pauli_mat = sqrt(2)^nqubits * pauli_basis_vec[symboltoint(pauli)+1]

    id = I(2^nqubits)

    U = cos(theta / 2) * id - 1.0im * sin(theta / 2) * pauli_mat

    return U
end


"""
    commutes(gate::PauliRotation, pstr::PauliString)

Check if a `PauliRotation` commutes with a `PauliString`.
"""
function commutes(gate::PauliRotation, pstr::PauliString)
    return commutes(gate, pstr.term)
end

"""
    commutes(gate::PauliRotation, pstr::Integer)

Check if a `PauliRotation` commutes with an integer Pauli string.
"""
function commutes(gate::PauliRotation, pstr::PauliStringType)
    # calculate the integer representation of the gate generator for faster computation
    generator_mask = symboltoint(typeof(pstr), gate.symbols, gate.qinds)
    return commutes(generator_mask, pstr)
end



# A parametrized Pauli rotation gate acting on the qubits `qinds` with the Pauli string `symbols`.
# The `term` is the integer representation of the Pauli string with the correct integer type for the total number of qubits.
# This allows for faster application of the gate.
# See `_tomaskedpaulirotation` for conversion from `PauliRotation`, which is the easiest way to construct a `MaskedPauliRotation`.
struct MaskedPauliRotation{T} <: ParametrizedGate where {T<:PauliStringType}
    symbols::Vector{Symbol}
    qinds::Vector{Int}
    generator_mask::T
end


# Union type for `PauliRotation` and `MaskedPauliRotation`, useful for functions which handle either agnostically.
PauliRotationUnion = Union{PauliRotation,MaskedPauliRotation}


# The fast `MaskedPauliRotation` version that we use in propagate()
function commutes(gate::MaskedPauliRotation, pstr::PauliStringType)
    return commutes(gate.generator_mask, pstr)
end


# Returns a circuit where `PauliRotation` gates are transformed to `MaskedPauliRotation` gates.
# This allows for significantly faster computation with the gate.
function _tomaskedpaulirotation(circ::Vector{G}, nqubits::Integer) where {G<:Gate}
    TT = getinttype(nqubits)
    return _tomaskedpaulirotation(circ, TT)
end


# Transforms a `PauliRotation` to a `MaskedPauliRotation` which carries the integer representation of the gate generator.
# This allows for significantly faster computation with the gate.
function _tomaskedpaulirotation(pauli_gate::PauliRotation, ::Type{TT}) where {TT<:PauliStringType}
    pstr_term = symboltoint(TT, pauli_gate.symbols, pauli_gate.qinds)
    return MaskedPauliRotation(pauli_gate.symbols, pauli_gate.qinds, pstr_term)
end


# Returns a circuit where `PauliRotation` gates are transformed to `MaskedPauliRotation` gates.
# This allows for significantly faster computation with the gate.
function _tomaskedpaulirotation(circ::Vector{G}, ::Type{TT}) where {G<:Gate,TT<:PauliStringType}
    masked_circ = copy(circ)
    for (ii, gate) in enumerate(masked_circ)
        if isa(gate, PauliRotation)
            masked_circ[ii] = _tomaskedpaulirotation(gate, TT)
        end
    end
    return masked_circ
end


# Transforms a `PauliRotation` to a `MaskedPauliRotation` which carries the integer representation of the gate generator.
# This allows for significantly faster computation with the gate.
function _tomaskedpaulirotation(pauli_gate::PauliRotation, nqubits::Integer)
    pstr_term = symboltoint(nqubits, pauli_gate.symbols, pauli_gate.qinds)
    return MaskedPauliRotation(pauli_gate.symbols, pauli_gate.qinds, pstr_term)
end


# Transforms a `FrozenGate` with a `PauliRotation` to a `MaskedPauliRotation` 
# which carries the integer representation of the gate generator.
function _tomaskedpaulirotation(frozen_gate::FrozenGate, nqubits::Integer)
    return FrozenGate(_tomaskedpaulirotation(frozen_gate.gate, nqubits), frozen_gate.theta)
end


# Return the `MaskedPauliRotation` gate as is.
function _tomaskedpaulirotation(masked_pauli_gate::MaskedPauliRotation, args...)
    return masked_pauli_gate
end




